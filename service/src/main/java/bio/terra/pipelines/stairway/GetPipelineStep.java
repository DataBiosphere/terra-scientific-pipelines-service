package bio.terra.pipelines.stairway;

import static java.lang.Boolean.TRUE;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
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
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put("updateComplete", Boolean.FALSE);

    Pipeline pipelineInfo = pipelinesService.getPipeline("imputation");
    workingMap.put("updateComplete", TRUE);
    workingMap.put(StairwayJobBuilder.JobMapKeys.RESPONSE.getKeyName(), pipelineInfo);
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    Boolean didUpdate = flightContext.getWorkingMap().get("updateComplete", Boolean.class);
    if (TRUE.equals(didUpdate)) {
      // If the update is complete, then we cannot undo it. This is a teeny tiny window
      // where the error occurs after the update, but before the success return. However,
      // the DebugInfo.lastStepFailure will always hit it.
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new RuntimeException("dismal failure"));
    }
    // Failed before update - perform undo
    return StepResult.getStepResultSuccess();
  }
}