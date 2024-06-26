package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.getBlobNameFromTerraWorkspaceStorageHttpUrl;
import static bio.terra.pipelines.common.utils.FileUtils.getFileNameFromFullPath;
import static bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService.READ_PERMISSION_STRING;
import static bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService.WRITE_PERMISSION_STRING;

import bio.terra.common.db.WriteTransaction;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutput;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlight;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
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
  private final SamService samService;
  private final PipelineRunsRepository pipelineRunsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;
  private final WorkspaceManagerService workspaceManagerService;

  private final ObjectMapper objectMapper = new ObjectMapper();

  @Autowired
  public PipelineRunsService(
      JobService jobService,
      SamService samService,
      PipelineRunsRepository pipelineRunsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      WorkspaceManagerService workspaceManagerService) {
    this.jobService = jobService;
    this.samService = samService;
    this.pipelineRunsRepository = pipelineRunsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
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

    PipelinesEnum pipelineName = PipelinesEnum.valueOf(pipeline.getName().toUpperCase());

    if (pipeline.getWorkspaceId() == null) {
      throw new InternalServerErrorException("%s workspaceId not defined".formatted(pipelineName));
    }

    // save the pipeline run to the database
    writePipelineRunToDb(
        jobId,
        userId,
        pipeline.getId(),
        pipeline.getWorkspaceId(),
        CommonPipelineRunStatusEnum.PREPARING,
        null,
        null,
        userProvidedInputs);

    // return a map of SAS urls and azcopy commands for the user to upload their input files
    Map<String, Map<String, String>> pipelineFileInputs =
        prepareFileInputs(pipeline, jobId, userProvidedInputs);

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
        pipeline.getPipelineInputDefinitions().stream()
            .filter(PipelineInputDefinition::getUserProvided)
            // note VCF is our only file type currently, and this method doesn't yet support
            // VCF_ARRAY
            .filter(p -> p.getType().equals(PipelineInputTypesEnum.VCF))
            .map(PipelineInputDefinition::getName)
            .toList();

    // generate a map where the key is the input name, and the value is a map containing the
    // write-only SAS url for the file and the full azcopy command to upload the file

    UUID workspaceId = pipeline.getWorkspaceId();
    UUID storageResourceId =
        workspaceManagerService.getWorkspaceStorageResourceId(
            workspaceId, samService.getTeaspoonsServiceAccountToken());
    String accessToken = samService.getTeaspoonsServiceAccountToken();

    Map<String, Map<String, String>> fileInputsMap = new HashMap<>();
    for (String fileInputName : fileInputNames) {
      fileInputsMap.put(
          fileInputName,
          constructSasUrlAndAzcopyCommand(
              jobId,
              workspaceId,
              storageResourceId,
              userProvidedInputs.get(fileInputName).toString(),
              accessToken));
    }

    return fileInputsMap;
  }

  private Map<String, String> constructSasUrlAndAzcopyCommand(
      UUID jobId,
      UUID workspaceId,
      UUID storageResourceId,
      String fileInputValue,
      String accessToken) {
    String userProvidedFileName = getFileNameFromFullPath(fileInputValue);
    String destinationBlobName = "user-input-files/%s/%s".formatted(jobId, userProvidedFileName);
    String sasUrl =
        workspaceManagerService.getSasUrlForBlob(
            workspaceId,
            storageResourceId,
            destinationBlobName,
            WRITE_PERMISSION_STRING,
            accessToken);
    Map<String, String> singleFileInputMap = new HashMap<>();
    singleFileInputMap.put("sasUrl", sasUrl);
    singleFileInputMap.put("azcopyCommand", "azcopy copy %s %s".formatted(fileInputValue, sasUrl));
    return singleFileInputMap;
  }

  /**
   * Create a new PipelineRun for a given pipeline, job, and user-provided inputs.
   *
   * <p>We encase the logic here in a transaction so that if the submission to Stairway fails, we do
   * not persist the entry in our own pipeline_runs table.
   *
   * <p>The Teaspoons database will auto-generate created and updated timestamps, so we do not need
   * to specify them when writing to the database. The generated timestamps will be included in the
   * returned PipelineRun object.
   */
  @WriteTransaction
  @SuppressWarnings("java:S1301") // allow switch statement with only one case
  public PipelineRun createPipelineRun(
      Pipeline pipeline,
      UUID jobId,
      String userId,
      String description,
      Map<String, Object> userProvidedInputs,
      String resultPath) {

    PipelinesEnum pipelineName = PipelinesEnum.valueOf(pipeline.getName().toUpperCase());

    if (pipeline.getWorkspaceId() == null) {
      throw new InternalServerErrorException("%s workspaceId not defined".formatted(pipelineName));
    }

    PipelineRun pipelineRun =
        writePipelineRunToDb(
            jobId,
            userId,
            pipeline.getId(),
            pipeline.getWorkspaceId(),
            CommonPipelineRunStatusEnum.SUBMITTED,
            description,
            resultPath,
            userProvidedInputs);

    logger.info("Creating new {} job for user {}", pipelineName, userId);

    Class<? extends Flight> flightClass;
    switch (pipelineName) {
      case IMPUTATION_BEAGLE:
        flightClass = RunImputationJobFlight.class;
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
            .addParameter(RunImputationJobFlightMapKeys.PIPELINE_ID, pipeline.getId())
            .addParameter(
                RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
                pipeline.getPipelineInputDefinitions())
            .addParameter(
                RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS, userProvidedInputs)
            .addParameter(
                RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID,
                pipeline.getWorkspaceId().toString())
            .addParameter(
                RunImputationJobFlightMapKeys.WDL_METHOD_NAME, pipeline.getWdlMethodName())
            .addParameter(JobMapKeys.RESULT_PATH.getKeyName(), resultPath);

    jobBuilder.submit();

    logger.info("Created {} pipelineRun with jobId {}", pipelineName, jobId);

    return pipelineRun;
  }

  @SuppressWarnings({"java:S107"}) // Disable "Methods should not have too many parameters"
  public PipelineRun writePipelineRunToDb(
      UUID jobUuid,
      String userId,
      Long pipelineId,
      UUID controlWorkspaceId,
      CommonPipelineRunStatusEnum status,
      String description,
      String resultUrl,
      Map<String, Object> pipelineInputs) {

    // write pipelineRun to database
    PipelineRun pipelineRun =
        new PipelineRun(
            jobUuid,
            userId,
            pipelineId,
            controlWorkspaceId,
            status.toString(),
            description,
            resultUrl);
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

  public PipelineRun getPipelineRun(UUID jobId, String userId) {
    return pipelineRunsRepository.findByJobIdAndUserId(jobId, userId).orElse(null);
  }

  /**
   * Mark a pipeline run as successful (is_success = True) in our database.
   *
   * <p>We expect this method to be called by the final step of a flight, at which point we assume
   * that the pipeline_run has completed successfully. Therefore, we do not do any checks on the
   * status column here. It is currently possible to mark an incomplete pipeline_run as is_success =
   * True using this method.
   */
  public PipelineRun markPipelineRunSuccessAndWriteOutputs(
      UUID jobId, String userId, Map<String, String> outputs) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);
    pipelineRun.setIsSuccess(true);
    pipelineRun.setOutput(pipelineRunOutputAsString(outputs));

    return pipelineRunsRepository.save(pipelineRun);
  }

  /**
   * Extract the pipeline outputs from a pipelineRun object, fetch SAS tokens for (currently all of)
   * them, and return an ApiPipelineRunOutput object with the formatted outputs.
   *
   * @param pipelineRun object from the pipelineRunsRepository
   * @return ApiPipelineRunOutput
   */
  public ApiPipelineRunOutput formatPipelineRunOutputs(PipelineRun pipelineRun) {
    Map<String, String> outputMap = pipelineRunOutputAsMap(pipelineRun.getOutput());

    UUID workspaceId = pipelineRun.getWorkspaceId();
    String accessToken = samService.getTeaspoonsServiceAccountToken();
    logger.info("Calling WSM to get storage container id for workspace: {}", workspaceId);
    UUID resourceId =
        workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, accessToken);

    // currently all outputs are paths that will need a SAS token
    outputMap.replaceAll(
        (k, v) ->
            workspaceManagerService.getSasUrlForBlob(
                workspaceId,
                resourceId,
                getBlobNameFromTerraWorkspaceStorageHttpUrl(v, workspaceId),
                READ_PERMISSION_STRING,
                accessToken));
    ApiPipelineRunOutput apiPipelineRunOutput = new ApiPipelineRunOutput();
    apiPipelineRunOutput.putAll(outputMap);
    return apiPipelineRunOutput;
  }

  public String pipelineRunOutputAsString(Map<String, String> outputMap) {
    try {
      return objectMapper.writeValueAsString(outputMap);
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error converting pipeline run output to string", e);
    }
  }

  public Map<String, String> pipelineRunOutputAsMap(String outputString) {
    try {
      return objectMapper.readValue(outputString, new TypeReference<>() {});
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error reading pipeline run outputs", e);
    }
  }
}
