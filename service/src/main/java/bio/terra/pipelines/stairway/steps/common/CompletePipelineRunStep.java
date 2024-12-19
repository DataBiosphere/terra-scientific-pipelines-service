package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to mark a pipeline run as a success and write the outputs and quota consumed to the database
 */
public class CompletePipelineRunStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final Logger logger = LoggerFactory.getLogger(CompletePipelineRunStep.class);

  public CompletePipelineRunStep(PipelineRunsService pipelineRunsService) {
    this.pipelineRunsService = pipelineRunsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID);

    UUID jobId = UUID.fromString(flightContext.getFlightId());
    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

    // validate and extract parameters from working map
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        workingMap, ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED);
    Map<String, String> outputsMap =
        workingMap.get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class);
    int quotaConsumed = workingMap.get(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, Integer.class);

    pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(
        jobId, userId, outputsMap, quotaConsumed);

    logger.info("Marked run {} as a success and wrote outputs to the db", jobId);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
