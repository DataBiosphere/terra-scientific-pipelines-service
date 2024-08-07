package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;
import static bio.terra.pipelines.common.utils.FileUtils.getBlobNameFromTerraWorkspaceStorageHttpUrl;
import static bio.terra.pipelines.common.utils.FileUtils.getStorageContainerUrlFromSasUrl;
import static java.util.Collections.emptyList;
import static org.springframework.data.domain.PageRequest.ofSize;

import bio.terra.common.db.WriteTransaction;
import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.pagination.CursorBasedPageable;
import bio.terra.pipelines.common.utils.pagination.FieldEqualsSpecification;
import bio.terra.pipelines.common.utils.pagination.PageResponse;
import bio.terra.pipelines.common.utils.pagination.PageSpecification;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputs;
import bio.terra.pipelines.stairway.imputation.RunImputationAzureJobFlight;
import bio.terra.pipelines.stairway.imputation.RunImputationAzureJobFlightMapKeys;
import bio.terra.stairway.Flight;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.HashMap;
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
  private final PipelinesService pipelinesService;
  private final SamService samService;
  private final PipelineRunsRepository pipelineRunsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;
  private final PipelineOutputsRepository pipelineOutputsRepository;
  private final WorkspaceManagerService workspaceManagerService;

  private final ObjectMapper objectMapper = new ObjectMapper();

  @Autowired
  public PipelineRunsService(
      JobService jobService,
      PipelinesService pipelinesService,
      SamService samService,
      PipelineRunsRepository pipelineRunsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      PipelineOutputsRepository pipelineOutputsRepository,
      WorkspaceManagerService workspaceManagerService) {
    this.jobService = jobService;
    this.pipelinesService = pipelinesService;
    this.samService = samService;
    this.pipelineRunsRepository = pipelineRunsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.pipelineOutputsRepository = pipelineOutputsRepository;
    this.workspaceManagerService = workspaceManagerService;
  }

  /**
   * Prepare a new PipelineRun for a given pipeline and user-provided inputs. The caller provides a
   * job uuid and any relevant pipeline inputs. Teaspoons writes the pipeline run to the database,
   * calls WSM to generate a write-only SAS url for each file input, and increments the pipeline
   * prepareRun counter metric.
   *
   * <p>Teaspoons returns a map of the pipeline inputs to the user, containing, for each file input,
   * the write-only SAS url for that file and the full azcopy command to upload the file using that
   * SAS url.
   *
   * @param pipeline the pipeline to run
   * @param jobId the job uuid
   * @param userId the user id
   * @param userProvidedInputs the user-provided inputs
   */
  @WriteTransaction
  public Map<String, Map<String, String>> preparePipelineRun(
      Pipeline pipeline, UUID jobId, String userId, Map<String, Object> userProvidedInputs) {

    PipelinesEnum pipelineName = pipeline.getName();

    if (pipeline.getWorkspaceId() == null) {
      throw new InternalServerErrorException("%s workspaceId not defined".formatted(pipelineName));
    }

    if (pipelineRunExistsWithJobId(jobId)) {
      throw new BadRequestException(
          "JobId %s already exists. If you submitted this job, you can use the getPipelineRunResult endpoint to see details for it."
              .formatted(jobId));
    }

    // return a map of SAS urls and azcopy commands for the user to upload their input files
    Map<String, Map<String, String>> pipelineFileInputs =
        prepareFileInputs(pipeline, jobId, userProvidedInputs);

    // get an arbitrary file input (since all file inputs are in the same workspace) and
    // extract its sasUrl to get the workspace storage URL
    String arbitraryInputFileSasUrl = pipelineFileInputs.values().iterator().next().get("sasUrl");
    String workspaceStorageContainerUrl =
        getStorageContainerUrlFromSasUrl(arbitraryInputFileSasUrl, pipeline.getWorkspaceId());

    // save the pipeline run to the database
    writeNewPipelineRunToDb(
        jobId,
        userId,
        pipeline.getId(),
        pipeline.getWorkspaceId(),
        pipeline.getWorkspaceProject(),
        pipeline.getWorkspaceName(),
        workspaceStorageContainerUrl,
        userProvidedInputs);

    // increment the prepare metric for this pipeline
    MetricsUtils.incrementPipelinePrepareRun(pipelineName);

    return pipelineFileInputs;
  }

  /**
   * Generate SAS urls and azcopy commands for each user-provided file input in the pipeline.
   *
   * <p>Each user-provided file input (assumed to be a path to a local file) is translated into a
   * write-only SAS url in a location in the pipeline workspace storage container, in a directory
   * defined by the jobId.
   *
   * <p>This SAS url along with the source file path provided by the user are used to generate an
   * azcopy command that the user can run to upload the file to the location in the pipeline
   * workspace storage container.
   */
  private Map<String, Map<String, String>> prepareFileInputs(
      Pipeline pipeline, UUID jobId, Map<String, Object> userProvidedInputs) {
    // get the list of files that the user needs to upload
    List<String> fileInputNames =
        pipelinesService.extractUserProvidedFileInputNames(pipeline.getPipelineInputDefinitions());

    // generate a map where the key is the input name, and the value is a map containing the
    // write-only SAS url for the file and the full azcopy command to upload the file

    UUID workspaceId = pipeline.getWorkspaceId();
    UUID storageResourceId =
        workspaceManagerService.getWorkspaceStorageResourceId(
            workspaceId, samService.getTeaspoonsServiceAccountToken());
    String accessToken = samService.getTeaspoonsServiceAccountToken();

    Map<String, Map<String, String>> fileInputsMap = new HashMap<>();
    for (String fileInputName : fileInputNames) {
      String fileInputValue = (String) userProvidedInputs.get(fileInputName);
      String sasUrl =
          retrieveWriteOnlySasUrl(
              jobId, workspaceId, storageResourceId, fileInputValue, accessToken);

      fileInputsMap.put(
          fileInputName,
          Map.of(
              "sasUrl",
              sasUrl,
              "azcopyCommand",
              "azcopy copy %s %s".formatted(fileInputValue, sasUrl)));
    }

    return fileInputsMap;
  }

  private String retrieveWriteOnlySasUrl(
      UUID jobId,
      UUID workspaceId,
      UUID storageResourceId,
      String fileInputValue,
      String accessToken) {
    String destinationBlobName =
        constructDestinationBlobNameForUserInputFile(jobId, fileInputValue);
    return workspaceManagerService.getWriteSasUrlForBlob(
        workspaceId, storageResourceId, destinationBlobName, accessToken);
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
  public PipelineRun startPipelineRun(
      Pipeline pipeline, UUID jobId, String userId, String description, String resultPath) {

    PipelinesEnum pipelineName = pipeline.getName();

    if (pipeline.getWorkspaceId() == null) {
      throw new InternalServerErrorException("%s workspaceId not defined".formatted(pipelineName));
    }

    PipelineRun pipelineRun = startPipelineRunInDb(jobId, userId, description, resultPath);

    Map<String, Object> userProvidedInputs = retrievePipelineInputs(pipelineRun);

    logger.info("Starting new {} job for user {}", pipelineName, userId);

    Class<? extends Flight> flightClass;
    switch (pipelineName) {
      case IMPUTATION_BEAGLE:
        flightClass = RunImputationAzureJobFlight.class;
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
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), pipelineName)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), userId)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description)
            .addParameter(RunImputationAzureJobFlightMapKeys.PIPELINE_ID, pipeline.getId())
            .addParameter(
                RunImputationAzureJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
                pipeline.getPipelineInputDefinitions())
            .addParameter(
                RunImputationAzureJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
                pipeline.getPipelineOutputDefinitions())
            .addParameter(
                RunImputationAzureJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
                userProvidedInputs)
            .addParameter(
                RunImputationAzureJobFlightMapKeys.CONTROL_WORKSPACE_ID,
                pipeline.getWorkspaceId().toString())
            .addParameter(
                RunImputationAzureJobFlightMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_URL,
                pipelineRun.getWorkspaceStorageContainerUrl())
            .addParameter(
                RunImputationAzureJobFlightMapKeys.WDL_METHOD_NAME, pipeline.getWdlMethodName())
            .addParameter(JobMapKeys.RESULT_PATH.getKeyName(), resultPath);

    jobBuilder.submit();

    logger.info("Started {} pipelineRun with jobId {}", pipelineName, jobId);

    return pipelineRun;
  }

  private Map<String, Object> retrievePipelineInputs(PipelineRun pipelineRun) {
    PipelineInput pipelineInput =
        pipelineInputsRepository
            .findById(pipelineRun.getId())
            .orElseThrow(
                () ->
                    new InternalServerErrorException(
                        "Pipeline inputs not found for jobId %s"
                            .formatted(pipelineRun.getJobId())));
    try {
      return objectMapper.readValue(pipelineInput.getInputs(), new TypeReference<>() {});
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error reading pipeline inputs", e);
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
      UUID controlWorkspaceId,
      String controlWorkspaceProject,
      String controlWorkspaceName,
      String workspaceStorageContainerUrl,
      Map<String, Object> pipelineInputs) {

    // write pipelineRun to database
    PipelineRun pipelineRun =
        new PipelineRun(
            jobUuid,
            userId,
            pipelineId,
            controlWorkspaceId,
            controlWorkspaceProject,
            controlWorkspaceName,
            workspaceStorageContainerUrl,
            CommonPipelineRunStatusEnum.PREPARING);
    PipelineRun createdPipelineRun = writePipelineRunToDbThrowsDuplicateException(pipelineRun);

    String pipelineInputsAsString;
    try {
      // do this to write the pipeline inputs without writing the class name
      pipelineInputsAsString = objectMapper.writeValueAsString(pipelineInputs);
    } catch (JsonProcessingException e) {
      // this should never happen
      throw new InternalServerErrorException("Internal error processing pipeline inputs", e);
    }

    // save related pipeline inputs
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(createdPipelineRun.getId());
    pipelineInput.setInputs(pipelineInputsAsString);
    pipelineInputsRepository.save(pipelineInput);

    return createdPipelineRun;
  }

  protected PipelineRun writePipelineRunToDbThrowsDuplicateException(PipelineRun pipelineRun)
      throws DuplicateObjectException {
    try {
      pipelineRunsRepository.save(pipelineRun);
      logger.info("pipelineRun saved for jobId: {}", pipelineRun.getJobId());
    } catch (DataIntegrityViolationException e) {
      if (e.getCause() instanceof ConstraintViolationException c
          && c.getConstraintName().contains("jobId_unique")) {
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
   * Mark a pipelineRun as RUNNING in our database and store the user-provided job description and
   * resultUrl.
   *
   * <p>We check that the pipelineRun already exists in our database and that the existing
   * pipelineRun has status PREPARING.
   *
   * @param jobId
   * @param userId
   * @param description
   * @param resultUrl
   * @return pipelineRun
   */
  public PipelineRun startPipelineRunInDb(
      UUID jobId, String userId, String description, String resultUrl) {
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
    pipelineRun.setDescription(description);
    pipelineRun.setResultUrl(resultUrl);

    return pipelineRunsRepository.save(pipelineRun);
  }

  /**
   * Mark a pipeline run as successful (is_success = True) in our database.
   *
   * <p>We expect this method to be called by the final step of a flight, at which point we assume
   * that the pipeline_run has completed successfully. Therefore, we do not do any checks on the
   * status column here. It is currently possible to mark an incomplete pipeline_run as is_success =
   * True using this method.
   */
  @WriteTransaction
  public PipelineRun markPipelineRunSuccessAndWriteOutputs(
      UUID jobId, String userId, Map<String, String> outputs) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(pipelineRun.getId());
    pipelineOutput.setOutputs(pipelineRunOutputsAsString(outputs));
    pipelineOutputsRepository.save(pipelineOutput);

    pipelineRun.setIsSuccess(true);

    return pipelineRunsRepository.save(pipelineRun);
  }

  // methods to interact with and format pipeline run outputs

  /**
   * Extract the pipeline outputs from a pipelineRun object, fetch SAS tokens for (currently all of)
   * them, and return an ApiPipelineRunOutputs object with the formatted outputs.
   *
   * @param pipelineRun object from the pipelineRunsRepository
   * @return ApiPipelineRunOutputs
   */
  public ApiPipelineRunOutputs formatPipelineRunOutputs(PipelineRun pipelineRun) {
    Map<String, String> outputsMap =
        pipelineRunOutputsAsMap(
            pipelineOutputsRepository.findPipelineOutputsByJobId(pipelineRun.getId()).getOutputs());

    UUID workspaceId = pipelineRun.getWorkspaceId();
    String accessToken = samService.getTeaspoonsServiceAccountToken();
    logger.info("Calling WSM to get storage container id for workspace: {}", workspaceId);
    UUID resourceId =
        workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, accessToken);

    // currently all outputs are paths that will need a SAS token
    outputsMap.replaceAll(
        (k, v) ->
            workspaceManagerService.getReadSasUrlForBlob(
                workspaceId,
                resourceId,
                getBlobNameFromTerraWorkspaceStorageHttpUrl(v, workspaceId),
                accessToken));
    ApiPipelineRunOutputs apiPipelineRunOutputs = new ApiPipelineRunOutputs();
    apiPipelineRunOutputs.putAll(outputsMap);
    return apiPipelineRunOutputs;
  }

  public String pipelineRunOutputsAsString(Map<String, String> outputsMap) {
    try {
      return objectMapper.writeValueAsString(outputsMap);
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error converting pipeline run outputs to string", e);
    }
  }

  public Map<String, String> pipelineRunOutputsAsMap(String outputsString) {
    try {
      return objectMapper.readValue(outputsString, new TypeReference<>() {});
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error reading pipeline run outputs", e);
    }
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
