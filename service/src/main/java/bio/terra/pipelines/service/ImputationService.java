package bio.terra.pipelines.service;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlight;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/** Service to encapsulate logic used to run an imputation pipeline */
@Service
public class ImputationService {
  private static final Logger logger = LoggerFactory.getLogger(ImputationService.class);
  private final PipelineRunsService pipelineRunsService;
  private final JobService jobService;

  @Autowired
  ImputationService(PipelineRunsService pipelineRunsService, JobService jobService) {
    this.pipelineRunsService = pipelineRunsService;
    this.jobService = jobService;
  }

  /**
   * Creates a new Imputation pipeline service job, using a Stairway flight, based on a user's
   * request. Returns jobId of new job (which is the same as the flightId) if flight submission is
   * successful.
   *
   * @param userId
   * @param imputationPipeline - a pipeline that handles imputation
   * @return String jobId
   *     <p>Note that the information in the requested job will grow over time, along with the
   *     following related classes:
   * @see PipelineRun
   */
  public UUID createImputationJob(
      UUID jobId,
      String userId,
      String description,
      Pipeline imputationPipeline,
      Map<String, Object> userProvidedPipelineInputs,
      String resultPath) {

    PipelinesEnum imputationPipelineName =
        PipelinesEnum.valueOf(imputationPipeline.getName().toUpperCase());
    logger.info("Create new {} job for user {}", imputationPipelineName, userId);

    JobBuilder jobBuilder =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(RunImputationJobFlight.class)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), userId)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description)
            .addParameter(RunImputationJobFlightMapKeys.PIPELINE_ID, imputationPipeline.getId())
            .addParameter(
                RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
                imputationPipeline.getPipelineInputDefinitions())
            .addParameter(
                RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
                userProvidedPipelineInputs)
            .addParameter(
                RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID,
                imputationPipeline.getWorkspaceId().toString())
            .addParameter(
                RunImputationJobFlightMapKeys.WDL_METHOD_NAME,
                imputationPipeline.getWdlMethodName())
            .addParameter(JobMapKeys.RESULT_PATH.getKeyName(), resultPath);

    return jobBuilder.submit();
  }
}
