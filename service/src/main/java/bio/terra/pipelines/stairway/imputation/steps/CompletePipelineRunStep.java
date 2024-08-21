package bio.terra.pipelines.stairway.imputation.steps;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID.getKeyName());

    UUID jobId = UUID.fromString(flightContext.getFlightId());
    String userId = inputParameters.get(JobMapKeys.USER_ID.getKeyName(), String.class);

    // validate and extract parameters from working map
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        workingMap, RunImputationJobFlightMapKeys.PIPELINE_RUN_OUTPUTS);
    Map<String, String> outputsMap =
        workingMap.get(RunImputationJobFlightMapKeys.PIPELINE_RUN_OUTPUTS, Map.class);

    pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(jobId, userId, outputsMap);

    logger.info("Marked run {} as a success and wrote outputs to the db", jobId);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
