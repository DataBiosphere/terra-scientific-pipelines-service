package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.constructFilePath;

import bio.terra.common.db.WriteTransaction;
import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.exception.ValidationException;
import bio.terra.common.iam.SamUser;
import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.app.configuration.external.GcsConfiguration;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineRunFilterSpecification;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.flights.datadelivery.DataDeliveryJobMapKeys;
import bio.terra.pipelines.stairway.flights.datadelivery.DeliverDataToGcsFlight;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.flights.imputation.v20251002.RunImputationGcpJobFlight;
import bio.terra.stairway.Flight;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.hibernate.exception.ConstraintViolationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.data.domain.Page;
import org.springframework.data.domain.PageRequest;
import org.springframework.data.domain.Sort;
import org.springframework.data.jpa.domain.Specification;
import org.springframework.stereotype.Service;

/** Service to encapsulate logic used to manage pipeline runs */
@Service
public class PipelineRunsService {
  private static final Logger logger = LoggerFactory.getLogger(PipelineRunsService.class);

  private final JobService jobService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final PipelineRunsRepository pipelineRunsRepository;
  private final IngressConfiguration ingressConfiguration;
  private final ToolConfigService toolConfigService;
  private final SamService samService;
  private final GcsConfiguration gcsConfiguration;

  public static final List<String> ALLOWED_SORT_PROPERTIES =
      List.of("created", "updated", "quotaConsumed");
  private final GcsService gcsService;

  @Autowired
  public PipelineRunsService(
      JobService jobService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      PipelineRunsRepository pipelineRunsRepository,
      IngressConfiguration ingressConfiguration,
      ToolConfigService toolConfigService,
      SamService samService,
      GcsConfiguration gcsConfiguration,
      GcsService gcsService) {
    this.jobService = jobService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.pipelineRunsRepository = pipelineRunsRepository;
    this.ingressConfiguration = ingressConfiguration;
    this.toolConfigService = toolConfigService;
    this.samService = samService;
    this.gcsConfiguration = gcsConfiguration;
    this.gcsService = gcsService;
  }

  /**
   * Prepare a new PipelineRun for a given pipeline and user-provided inputs. The caller provides a
   * job uuid and any relevant pipeline inputs. Teaspoons writes the pipeline run to the database
   * and increments the pipeline prepareRun counter metric.
   *
   * <p>If the request contains local file inputs, Teaspoons returns a map of the pipeline file
   * inputs to the user, containing a signed URL string and a curl command they can use to upload
   * each input.
   *
   * <p>If the request contains cloud-based file inputs, Teaspoons checks that the user and the
   * Teaspoons SA has access to the bucket containing the file inputs.
   *
   * @param pipeline the pipeline to run
   * @param jobId the job uuid
   * @param authedUser the sam authedUser object for the user
   * @param userProvidedInputs the user-provided inputs
   * @return if local inputs to upload, a map of pipeline file inputs containing signed URLs and
   *     curl commands for the user to upload their files. if cloud inputs, return nothing.
   */
  @WriteTransaction
  public Map<String, Map<String, String>> preparePipelineRunV2(
      Pipeline pipeline,
      UUID jobId,
      SamUser authedUser,
      Map<String, Object> userProvidedInputs,
      String description,
      Boolean useResumableUploads) {

    validatePipelineWorkspaceSetup(pipeline);

    if (pipelineRunExistsWithJobId(jobId)) {
      throw new BadRequestException(
          "JobId %s already exists. If you submitted this job, you can use the getPipelineRunResult endpoint to see details for it."
              .formatted(jobId));
    }

    Map<String, Map<String, String>> pipelineFileInputSignedUrls;
    if (pipelineInputsOutputsService.userProvidedInputsAreGcsCloud(pipeline, userProvidedInputs)) {
      // do user and service access checks
      logger.info(
          "Found cloud inputs for jobId {}, no signed URLs needed; checking read access to input files",
          jobId);
      pipelineInputsOutputsService.validateUserAndServiceReadAccessToCloudInputs(
          pipeline, userProvidedInputs, authedUser);
      pipelineFileInputSignedUrls = null;
    } else {
      // return a map of signed PUT urls and curl commands for the user to upload their input files
      pipelineFileInputSignedUrls =
          pipelineInputsOutputsService.prepareLocalFileInputs(
              pipeline, jobId, userProvidedInputs, useResumableUploads);
    }

    // add default values to any optional inputs not specified by the user
    Map<String, Object> userProvidedInputsWithDefaults =
        pipelineInputsOutputsService.populateDefaultValuesForMissingOptionalUserInputs(
            pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    // save the pipeline run to the database
    writeNewPipelineRunToDb(
        jobId,
        authedUser.getSubjectId(),
        pipeline.getId(),
        pipeline.getToolVersion(),
        pipeline.getWorkspaceBillingProject(),
        pipeline.getWorkspaceName(),
        pipeline.getWorkspaceStorageContainerName(),
        pipeline.getWorkspaceGoogleProject(),
        userProvidedInputsWithDefaults,
        description);

    // increment the prepare metric for this pipeline
    MetricsUtils.incrementPipelinePrepareRun(pipeline.getName());

    return pipelineFileInputSignedUrls;
  }

  /**
   * Prepare a new PipelineRun for a given pipeline and user-provided inputs. The caller provides a
   * job uuid and any relevant pipeline inputs. Teaspoons writes the pipeline run to the database
   * and increments the pipeline prepareRun counter metric.
   *
   * <p>Teaspoons returns a map of the pipeline file inputs to the user, containing a signed URL
   * string and a curl command they can use to upload each input.
   *
   * @param pipeline the pipeline to run
   * @param jobId the job uuid
   * @param userId the user id
   * @param userProvidedInputs the user-provided inputs
   * @return a map of pipeline file inputs containing signed URLs and curl commands for the user to
   *     upload their files
   * @deprecated
   */
  @Deprecated(since = "2.2.0")
  @WriteTransaction
  public Map<String, Map<String, String>> preparePipelineRun(
      Pipeline pipeline,
      UUID jobId,
      String userId,
      Map<String, Object> userProvidedInputs,
      String description,
      Boolean useResumableUploads) {

    validatePipelineWorkspaceSetup(pipeline);

    if (pipelineRunExistsWithJobId(jobId)) {
      throw new BadRequestException(
          "JobId %s already exists. If you submitted this job, you can use the getPipelineRunResult endpoint to see details for it."
              .formatted(jobId));
    }

    // create signed PUT urls and curl commands for the user to upload their input files
    Map<String, Map<String, String>> pipelineFileInputSignedUrls =
        pipelineInputsOutputsService.prepareLocalFileInputs(
            pipeline, jobId, userProvidedInputs, useResumableUploads);

    // add default values to any optional inputs not specified by the user
    userProvidedInputs =
        pipelineInputsOutputsService.populateDefaultValuesForMissingOptionalUserInputs(
            pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    // save the pipeline run to the database
    writeNewPipelineRunToDb(
        jobId,
        userId,
        pipeline.getId(),
        pipeline.getToolVersion(),
        pipeline.getWorkspaceBillingProject(),
        pipeline.getWorkspaceName(),
        pipeline.getWorkspaceStorageContainerName(),
        pipeline.getWorkspaceGoogleProject(),
        userProvidedInputs,
        description);

    // increment the prepare metric for this pipeline
    MetricsUtils.incrementPipelinePrepareRun(pipeline.getName());

    return pipelineFileInputSignedUrls;
  }

  /**
   * Start a PipelineRun that exists in the database (via preparePipelineRun).
   *
   * <p>We encase the logic here in a transaction so that if the submission to Stairway fails, we do
   * not update the status in our own pipeline_runs table.
   *
   * <p>The Teaspoons database will auto-generate updated timestamps.
   */
  @WriteTransaction
  @SuppressWarnings("java:S1301") // allow switch statement with only one case
  public PipelineRun startPipelineRun(Pipeline pipeline, UUID jobId, String userId) {

    PipelinesEnum pipelineName = pipeline.getName();

    validatePipelineWorkspaceSetup(pipeline);

    PipelineRun pipelineRun = startPipelineRunInDb(jobId, userId);

    Map<String, Object> userProvidedInputs =
        pipelineInputsOutputsService.retrieveUserProvidedInputs(pipelineRun);

    logger.info("Starting new {} job for user {}", pipelineName, userId);

    Class<? extends Flight> flightClass;
    switch (pipelineName) {
      case ARRAY_IMPUTATION:
        flightClass = RunImputationGcpJobFlight.class; // v20251002
        break;
      default:
        throw new InternalServerErrorException(
            "Pipeline %s not supported by PipelineRunsService".formatted(pipelineName));
    }

    JobBuilder jobBuilder =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(flightClass)
            .addParameter(JobMapKeys.PIPELINE_NAME, pipelineName)
            .addParameter(JobMapKeys.PIPELINE_VERSION, pipeline.getVersion())
            .addParameter(JobMapKeys.USER_ID, userId)
            .addParameter(JobMapKeys.DESCRIPTION, pipelineRun.getDescription())
            .addParameter(JobMapKeys.PIPELINE_ID, pipeline.getId())
            .addParameter(JobMapKeys.DOMAIN_NAME, ingressConfiguration.getDomainName())
            .addParameter(JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK, true)
            .addParameter(JobMapKeys.DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK, true)
            .addParameter(JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK, true)
            .addParameter(ImputationJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS, userProvidedInputs)
            .addParameter(
                ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
                pipeline.getWorkspaceBillingProject())
            .addParameter(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, pipeline.getWorkspaceName())
            .addParameter(
                ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
                pipelineRun.getWorkspaceStorageContainerName())
            .addParameter(
                ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL,
                "gs://") // this is the GCP storage url protocol
            .addParameter(
                ImputationJobMapKeys.PIPELINE_TOOL_CONFIG,
                toolConfigService.getPipelineMainToolConfig(pipeline))
            .addParameter(
                ImputationJobMapKeys.QUOTA_TOOL_CONFIG,
                toolConfigService.getQuotaConsumedToolConfig(pipeline))
            .addParameter(
                ImputationJobMapKeys.INPUT_QC_TOOL_CONFIG,
                toolConfigService.getInputQcToolConfig(pipeline));

    jobBuilder.submit();

    logger.info("Started {} pipelineRun with jobId {}", pipelineName, jobId);

    return pipelineRun;
  }

  /** Validate that the pipeline object has workspace fields defined. */
  private static void validatePipelineWorkspaceSetup(Pipeline pipeline) {
    if (pipeline.getWorkspaceBillingProject() == null
        || pipeline.getWorkspaceName() == null
        || pipeline.getWorkspaceStorageContainerName() == null) {
      throw new InternalServerErrorException(
          "%s workspace not defined".formatted(pipeline.getName()));
    }
  }

  // methods to write and update PipelineRuns in the database

  /**
   * Write a new pipelineRun to the database, including the pipeline inputs. Its status will be set
   * to PREPARING.
   *
   * <p>The Teaspoons database will auto-generate created and updated timestamps, so we do not need
   * to specify them when writing to the database. The generated timestamps will be included in the
   * returned PipelineRun object.
   */
  @SuppressWarnings({"java:S107"}) // Disable "Methods should not have too many parameters"
  public PipelineRun writeNewPipelineRunToDb(
      UUID jobUuid,
      String userId,
      Long pipelineId,
      String toolVersion,
      String controlWorkspaceProject,
      String controlWorkspaceName,
      String controlWorkspaceStorageContainerUrl,
      String controlWorkspaceGoogleProject,
      Map<String, Object> pipelineInputs,
      String description) {

    // write pipelineRun to database
    PipelineRun pipelineRun =
        new PipelineRun(
            jobUuid,
            userId,
            pipelineId,
            toolVersion,
            controlWorkspaceProject,
            controlWorkspaceName,
            controlWorkspaceStorageContainerUrl,
            controlWorkspaceGoogleProject,
            CommonPipelineRunStatusEnum.PREPARING,
            description);
    PipelineRun createdPipelineRun = writePipelineRunToDbThrowsDuplicateException(pipelineRun);

    pipelineInputsOutputsService.savePipelineInputs(createdPipelineRun.getId(), pipelineInputs);

    return createdPipelineRun;
  }

  protected PipelineRun writePipelineRunToDbThrowsDuplicateException(PipelineRun pipelineRun)
      throws DuplicateObjectException {
    try {
      pipelineRunsRepository.save(pipelineRun);
      logger.info("pipelineRun saved for jobId: {}", pipelineRun.getJobId());
    } catch (DataIntegrityViolationException e) {
      if (e.getCause() instanceof ConstraintViolationException c
          && c.getConstraintName().contains("pipeline_runs_jobId_unique")) {
        throw new DuplicateObjectException(
            String.format("Duplicate jobId %s found", pipelineRun.getJobId()));
      }
      throw e;
    }

    return pipelineRun;
  }

  private boolean pipelineRunExistsWithJobId(UUID jobId) {
    return pipelineRunsRepository.existsByJobId(jobId);
  }

  public PipelineRun getPipelineRun(UUID jobId, String userId) {
    return pipelineRunsRepository.findByJobIdAndUserId(jobId, userId).orElse(null);
  }

  public UUID submitDataDeliveryFlight(
      PipelineRun pipelineRun, UUID deliveryJobId, String destinationPath, SamUser authedUser) {
    GcsFile fullPathWithJobId =
        new GcsFile(constructFilePath(destinationPath, pipelineRun.getJobId().toString()));

    validateUserAndServiceWriteAccessToDestinationPath(fullPathWithJobId, authedUser);

    JobBuilder jobBuilder =
        jobService
            .newJob()
            .jobId(deliveryJobId)
            .flightClass(DeliverDataToGcsFlight.class)
            .addParameter(JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK, false)
            .addParameter(JobMapKeys.DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK, false)
            .addParameter(JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK, false)
            .addParameter(JobMapKeys.USER_ID, authedUser.getSubjectId())
            .addParameter(JobMapKeys.DOMAIN_NAME, ingressConfiguration.getDomainName())
            .addParameter(JobMapKeys.PIPELINE_NAME, pipelineRun.getPipeline().getName())
            .addParameter(JobMapKeys.PIPELINE_ID, pipelineRun.getPipeline().getId())
            .addParameter(
                JobMapKeys.DESCRIPTION, "Data delivery for pipeline run " + pipelineRun.getId())
            .addParameter(DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, fullPathWithJobId)
            .addParameter(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, pipelineRun.getJobId());

    jobBuilder.submit();

    logger.info(
        "Started data delivery flight {} for pipeline run {} to destination {}",
        deliveryJobId,
        pipelineRun.getId(),
        destinationPath);

    return deliveryJobId;
  }

  /**
   * Mark a pipelineRun as RUNNING in our database.
   *
   * <p>We check that the pipelineRun already exists in our database and that the existing
   * pipelineRun has status PREPARING.
   *
   * @param jobId
   * @param userId
   * @return pipelineRun
   */
  public PipelineRun startPipelineRunInDb(UUID jobId, String userId) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);
    if (pipelineRun == null) {
      throw new BadRequestException(
          "JobId %s not found. You must prepare a pipeline run before starting it."
              .formatted(jobId));
    }
    // only allow starting a pipeline run if it is in the PREPARING state
    if (!pipelineRun.getStatus().equals(CommonPipelineRunStatusEnum.PREPARING)) {
      throw new BadRequestException(
          "JobId %s is not in the PREPARING state. Cannot start pipeline run.".formatted(jobId));
    }
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.RUNNING);

    return pipelineRunsRepository.save(pipelineRun);
  }

  /** Set the raw quota consumed for a pipeline run in our database. */
  public void setPipelineRunRawQuotaConsumed(UUID jobId, String userId, int rawQuotaConsumed) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);
    pipelineRun.setRawQuotaConsumed(rawQuotaConsumed);
    pipelineRunsRepository.save(pipelineRun);
  }

  public void validateUserAndServiceWriteAccessToDestinationPath(
      GcsFile destinationPath, SamUser authedUser) {

    boolean userHasBucketWriteAccess =
        gcsService.userHasBucketWriteAccess(
            destinationPath.getBucketName(),
            samService.getUserPetServiceAccountTokenReadOnly(authedUser).getToken());

    if (!userHasBucketWriteAccess) {
      String userProxyGroup = samService.getProxyGroupForUser(authedUser);
      throw new ValidationException(
          "User %s does not have necessary permissions to write to destination bucket %s, or the bucket does not exist. Please ensure the user's proxy group %s has write access to the destination bucket."
              .formatted(authedUser, destinationPath.getBucketName(), userProxyGroup));
    }

    boolean serviceHasBucketWriteAccess =
        gcsService.serviceHasBucketWriteAccess(destinationPath.getBucketName());

    if (!serviceHasBucketWriteAccess) {
      throw new ValidationException(
          "Service does not have necessary permissions to write to destination bucket %s, or the bucket does not exist. Please ensure that %s has write access to the destination bucket."
              .formatted(
                  destinationPath.getBucketName(),
                  gcsConfiguration.serviceAccountGroupForCloudIntegration()));
    }
  }

  /**
   * Mark a pipeline run as successful (status = SUCCEEDED) in our database and store the effective
   * quota consumed by the job.
   *
   * <p>We expect this method to be called by the final step of a flight, at which point we assume
   * that the pipeline_run has completed successfully. Therefore, we do not do any checks on the
   * status column here. It is currently possible to mark an incomplete pipeline_run as is_success =
   * True using this method.
   */
  @WriteTransaction
  public PipelineRun markPipelineRunSuccessAndWriteOutputs(
      UUID jobId, String userId, Map<String, String> outputs, int quotaConsumed) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);

    pipelineInputsOutputsService.savePipelineOutputs(pipelineRun.getId(), outputs);

    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRun.setQuotaConsumed(quotaConsumed);

    return pipelineRunsRepository.save(pipelineRun);
  }

  /**
   * Mark a pipeline run as failed (status = FAILED) in our database.
   *
   * <p>We expect this method to be called by the undoStep method of the first step in a flight, so
   * that it is executed when the flight has failed.
   */
  public void markPipelineRunFailed(UUID jobId, String userId) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.FAILED);
    pipelineRunsRepository.save(pipelineRun);
  }

  /**
   * Extract a paginated list of Pipeline Run records from the database with filtering support
   *
   * @param pageNumber - the page number to retrieve
   * @param pageSize - how many records to return
   * @param sortProperty - which property to sort on
   * @param sortDirection - which direction to sort
   * @param userId - caller's user id
   * @param filters - map of field names to filter values
   * @return - a Page containing the list of records in the current page
   */
  public Page<PipelineRun> findPipelineRunsPaginated(
      int pageNumber,
      int pageSize,
      String sortProperty,
      String sortDirection,
      String userId,
      Map<String, String> filters) {
    // Validate and/or set defaults for the sort property and direction
    String validatedSortProperty = validateSortProperty(sortProperty);
    Sort.Direction validatedSortDirection = validateSortDirection(sortDirection);

    // Build the specification, which also validates the filters
    Specification<PipelineRun> filterSpec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, userId);

    PageRequest pageRequest =
        PageRequest.of(pageNumber, pageSize, validatedSortDirection, validatedSortProperty);

    return pipelineRunsRepository.findAll(filterSpec, pageRequest);
  }

  /**
   * Extract a paginated list of Pipeline Run records from the database
   *
   * @param pageNumber - the page number to retrieve
   * @param pageSize - how many records to return
   * @param sortProperty - which property to sort on
   * @param sortDirection - which direction to sort
   * @param userId - caller's user id
   * @return - a Page containing the list of records in the current page
   */
  public Page<PipelineRun> findPipelineRunsPaginated(
      int pageNumber, int pageSize, String sortProperty, String sortDirection, String userId) {
    return findPipelineRunsPaginated(
        pageNumber, pageSize, sortProperty, sortDirection, userId, Collections.emptyMap());
  }

  // Returns the total count of pipeline runs for a user
  public long getPipelineRunCount(String userId) {
    return pipelineRunsRepository.countByUserId(userId);
  }

  // Returns the total count of pipeline runs for a user with filters applied
  public long getFilteredPipelineRunCount(String userId, Map<String, String> filters) {
    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, userId);

    return pipelineRunsRepository.count(spec);
  }

  private Sort.Direction validateSortDirection(String sortDirection) {
    if (sortDirection == null) {
      return Sort.Direction.DESC;
    }
    try {
      return Sort.Direction.valueOf(sortDirection.toUpperCase());
    } catch (IllegalArgumentException e) {
      throw new BadRequestException(
          "Invalid sort direction: %s. Valid values are ASC or DESC".formatted(sortDirection));
    }
  }

  private String validateSortProperty(String sortProperty) {
    // default to sorting by created
    if (sortProperty == null) {
      return "created";
    }

    if (!ALLOWED_SORT_PROPERTIES.contains(sortProperty)) {
      throw new BadRequestException(
          "Invalid sort property: %s. Valid sort properties are: %s"
              .formatted(sortProperty, String.join(", ", ALLOWED_SORT_PROPERTIES)));
    }

    return sortProperty;
  }
}
