package bio.terra.pipelines.dependencies.stairway;

import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import org.springframework.http.HttpStatus;

public class StairwayJobServiceTestStep implements Step {

  public StairwayJobServiceTestStep() {}

  @Override
  public StepResult doStep(FlightContext context) {
    // Pull description from input parameters
    String description = context.getInputParameters().get("description", String.class);

    // Configure the results to be the description from inputs
    context.getWorkingMap().put(StairwayJobMapKeys.RESPONSE.getKeyName(), description);
    context
        .getWorkingMap()
        .put(StairwayJobMapKeys.STATUS_CODE.getKeyName(), HttpStatus.I_AM_A_TEAPOT);
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    return StepResult.getStepResultSuccess();
  }
}
