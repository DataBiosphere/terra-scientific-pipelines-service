package bio.terra.pipelines.stairway;

import static java.lang.Boolean.TRUE;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GetPipelineStep implements Step {

  private final PipelinesService pipelinesService;
  private final Logger logger = LoggerFactory.getLogger(GetPipelineStep.class);

  public GetPipelineStep(PipelinesService pipelinesService) {
    this.pipelinesService = pipelinesService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, GetPipelineFlightMapKeys.PIPELINE_ID);

    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(GetPipelineFlightMapKeys.LOOKUP_COMPLETE, Boolean.FALSE);

    Pipeline pipelineInfo =
        pipelinesService.getPipeline(
            inputParameters.get(GetPipelineFlightMapKeys.PIPELINE_ID, String.class));
    logger.info("Retrieved pipeline info: {}", pipelineInfo);

    workingMap.put(GetPipelineFlightMapKeys.LOOKUP_COMPLETE, TRUE);
    workingMap.put(StairwayJobMapKeys.RESPONSE.getKeyName(), pipelineInfo);

    FlightUtils.validateRequiredEntries(
        workingMap,
        GetPipelineFlightMapKeys.LOOKUP_COMPLETE,
        StairwayJobMapKeys.RESPONSE.getKeyName());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    Boolean didLookup =
        flightContext.getWorkingMap().get(GetPipelineFlightMapKeys.LOOKUP_COMPLETE, Boolean.class);
    if (TRUE.equals(didLookup)) {
      // If the action performed by doStep (database lookup in this toy example) is complete, then
      // we cannot undo it. (This is meaningless for a lookup but will make sense for real steps.)
      // This is a teeny tiny window
      // where the error occurs after the action, but before the success return. However,
      // the DebugInfo.lastStepFailure will always hit it.
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new RuntimeException("dismal failure"));
    }
    // Failed before update - perform undo - nothing to do here
    return StepResult.getStepResultSuccess();
  }
}
