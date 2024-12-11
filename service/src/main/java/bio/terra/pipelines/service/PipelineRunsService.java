package bio.terra.pipelines.service;

import static java.util.Collections.emptyList;
import static org.springframework.data.domain.PageRequest.ofSize;

import bio.terra.common.db.WriteTransaction;
import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.pagination.CursorBasedPageable;
import bio.terra.pipelines.common.utils.pagination.FieldEqualsSpecification;
import bio.terra.pipelines.common.utils.pagination.PageResponse;
import bio.terra.pipelines.common.utils.pagination.PageSpecification;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.flights.imputation.RunImputationGcpJobFlight;
import bio.terra.stairway.Flight;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.hibernate.exception.ConstraintViolationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.stereotype.Service;

/** Service to encapsulate logic used to manage pipeline runs */
@Service
public class PipelineRunsService {
  private static final Logger logger = LoggerFactory.getLogger(PipelineRunsService.class);

  private final JobService jobService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final PipelineRunsRepository pipelineRunsRepository;
  private final IngressConfiguration ingressConfiguration;

  @Autowired
  public PipelineRunsService(
      JobService jobService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      PipelineRunsRepository pipelineRunsRepository,
      IngressConfiguration ingressConfiguration) {
    this.jobService = jobService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.pipelineRunsRepository = pipelineRunsRepository;
    this.ingressConfiguration = ingressConfiguration;
  }

  /**
   * Prepare a new PipelineRun for a given pipeline and user-provided inputs. The caller provides a
   * job uuid and any relevant pipeline inputs. Teaspoons writes the pipeline run to the database
   * and increments the pipeline prepareRun counter metric.
   *
   * <p>Teaspoons returns a map of the pipeline inputs to the user, containing, for each file input,
   * path that they provided. Note that in future this will change depending on the GCP data
   * strategy we choose.
   *
   * @param pipeline the pipeline to run
   * @param jobId the job uuid
   * @param userId the user id
   * @param userProvidedInputs the user-provided inputs
   */
  @WriteTransaction
  public Map<String, Map<String, String>> preparePipelineRun(
      Pipeline pipeline,
      UUID jobId,
      String userId,
      Map<String, Object> userProvidedInputs,
      String description) {

    PipelinesEnum pipelineName = pipeline.getName();

    if (pipeline.getWorkspaceBillingProject() == null
        || pipeline.getWorkspaceName() == null
        || pipeline.getWorkspaceStorageContainerName() == null) {
      throw new InternalServerErrorException("%s workspace not defined".formatted(pipelineName));
    }

    if (pipelineRunExistsWithJobId(jobId)) {
      throw new BadRequestException(
          "JobId %s already exists. If you submitted this job, you can use the getPipelineRunResult endpoint to see details for it."
              .formatted(jobId));
    }

    // return a map of signed PUT urls and curl commands for the user to upload their input files
    Map<String, Map<String, String>> pipelineFileInputs =
        pipelineInputsOutputsService.prepareFileInputs(pipeline, jobId, userProvidedInputs);

    // save the pipeline run to the database
    writeNewPipelineRunToDb(
        jobId,
        userId,
        pipeline.getId(),
        pipeline.getWdlMethodVersion(),
        pipeline.getWorkspaceBillingProject(),
        pipeline.getWorkspaceName(),
        pipeline.getWorkspaceStorageContainerName(),
        pipeline.getWorkspaceGoogleProject(),
        userProvidedInputs,
        description);

    // increment the prepare metric for this pipeline
    MetricsUtils.incrementPipelinePrepareRun(pipelineName);

    return pipelineFileInputs;
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

    if (pipeline.getWorkspaceBillingProject() == null
        || pipeline.getWorkspaceName() == null
        || pipeline.getWorkspaceStorageContainerName() == null) {
      throw new InternalServerErrorException("%s workspace not defined".formatted(pipelineName));
    }

    PipelineRun pipelineRun = startPipelineRunInDb(jobId, userId);

    Map<String, Object> userProvidedInputs =
        pipelineInputsOutputsService.retrievePipelineInputs(pipelineRun);

    logger.info("Starting new {} job for user {}", pipelineName, userId);

    Class<? extends Flight> flightClass;
    switch (pipelineName) {
      case ARRAY_IMPUTATION:
        flightClass = RunImputationGcpJobFlight.class;
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
            .addParameter(JobMapKeys.USER_ID, userId)
            .addParameter(JobMapKeys.DESCRIPTION, pipelineRun.getDescription())
            .addParameter(JobMapKeys.PIPELINE_ID, pipeline.getId())
            .addParameter(JobMapKeys.DOMAIN_NAME, ingressConfiguration.getDomainName())
            .addParameter(JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK, true)
            .addParameter(JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK, true)
            .addParameter(
                ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
                pipeline.getPipelineInputDefinitions())
            .addParameter(
                ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
                pipeline.getPipelineOutputDefinitions())
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
            .addParameter(ImputationJobMapKeys.WDL_METHOD_NAME, pipeline.getWdlMethodName())
            .addParameter(ImputationJobMapKeys.WDL_METHOD_VERSION, pipeline.getWdlMethodVersion());

    jobBuilder.submit();

    logger.info("Started {} pipelineRun with jobId {}", pipelineName, jobId);

    return pipelineRun;
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
      String wdlMethodVersion,
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
            wdlMethodVersion,
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

  /**
   * Mark a pipeline run as successful (status = SUCCEEDED) in our database and store the quota
   * consumed by the job.
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
   * Extract a paginated list of Pipeline Run records from the database
   *
   * @param limit - how many records to return
   * @param pageToken - encoded token representing where to start the cursor based pagination from
   * @param userId - caller's user id
   * @return - a PageResponse containing the list of records in the current page and the page tokens
   *     for the next and previous page if applicable
   */
  public PageResponse<List<PipelineRun>> findPipelineRunsPaginated(
      int limit, String pageToken, String userId) {

    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(limit, pageToken, null);
    PageSpecification<PipelineRun> pageSpecification =
        new PageSpecification<>("id", cursorBasedPageable);
    FieldEqualsSpecification<PipelineRun> userIdSpecification =
        new FieldEqualsSpecification<>("userId", userId);

    var postSlice =
        pipelineRunsRepository.findAll(
            userIdSpecification.and(pageSpecification), ofSize(cursorBasedPageable.getSize()));
    if (!postSlice.hasContent()) return new PageResponse<>(emptyList(), null, null);

    var pipelineRuns = postSlice.getContent();
    return new PageResponse<>(
        pipelineRuns,
        CursorBasedPageable.getEncodedCursor(
            pipelineRuns.get(0).getId().toString(),
            pipelineRunsRepository.existsByIdGreaterThan(pipelineRuns.get(0).getId())),
        CursorBasedPageable.getEncodedCursor(
            pipelineRuns.get(pipelineRuns.size() - 1).getId().toString(), postSlice.hasNext()));
  }
}
