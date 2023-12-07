package bio.terra.pipelines.dependencies.stairway;

import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import org.springframework.http.HttpStatus;

public class StairwayJobServiceTestStep implements Step {
  private final String description;

  public StairwayJobServiceTestStep(String description) {
    this.description = description;
  }

  @Override
  public StepResult doStep(FlightContext context) {
    // Configure the results
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
