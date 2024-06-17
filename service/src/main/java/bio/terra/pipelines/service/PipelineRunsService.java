package bio.terra.pipelines.service;

import static bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService.READ_PERMISSION_STRING;

import bio.terra.common.db.WriteTransaction;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
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
import com.fasterxml.jackson.databind.ObjectMapper;
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
   * Create a new PipelineRun for a given pipeline, job, and user-provided inputs.
   *
   * <p>We encase the logic here in a transaction so that if the submission to Stairway fails, we do
   * not persist the entry in our own pipeline_runs table.
   *
   * <p>The TSPS database will auto-generate created and updated timestamps, so we do not need to
   * specify them when writing to the database. The generated timestamps will be included in the
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
    pipelineRun.setOutput(outputs);

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
    Map<String, String> outputMap = pipelineRun.getOutput();

    // currently all outputs are paths that will need a SAS token
    outputMap.replaceAll(
        (k, v) ->
            workspaceManagerService.getSasTokenForFile(
                pipelineRun.getWorkspaceId(),
                v,
                READ_PERMISSION_STRING,
                samService.getTspsServiceAccountToken()));
    ApiPipelineRunOutput apiPipelineRunOutput = new ApiPipelineRunOutput();
    apiPipelineRunOutput.putAll(outputMap);
    return apiPipelineRunOutput;
  }
}
