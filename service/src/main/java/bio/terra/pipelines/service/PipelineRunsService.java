package bio.terra.pipelines.service;

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
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
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
  private final PipelineRunsRepository pipelineRunsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;

  private final ObjectMapper objectMapper = new ObjectMapper();

  @Autowired
  public PipelineRunsService(
      JobService jobService,
      PipelineRunsRepository pipelineRunsRepository,
      PipelineInputsRepository pipelineInputsRepository) {
    this.jobService = jobService;
    this.pipelineRunsRepository = pipelineRunsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
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
      CommonPipelineRunStatusEnum status,
      String description,
      String resultUrl,
      Map<String, Object> pipelineInputs) {

    // write pipelineRun to database
    PipelineRun pipelineRun =
        new PipelineRun(jobUuid, userId, pipelineId, status.toString(), description, resultUrl);
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
  public PipelineRun markPipelineRunSuccess(UUID jobId, String userId) {
    PipelineRun pipelineRun = getPipelineRun(jobId, userId);
    pipelineRun.setIsSuccess(true);
    return pipelineRunsRepository.save(pipelineRun);
  }
}