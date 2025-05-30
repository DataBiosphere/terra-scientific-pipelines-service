package bio.terra.pipelines.stairway.steps.imputation.azure;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.stairway.*;

/** This step queries the Leonardo status endpoint to check if it is healthy */
public class CheckLeonardoHealthStep implements Step {
  private final LeonardoService leonardoService;

  public CheckLeonardoHealthStep(LeonardoService leonardoService) {
    this.leonardoService = leonardoService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    HealthCheck.Result healthResult = leonardoService.checkHealth();
    if (!healthResult.isOk()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_RETRY,
          new LeonardoServiceApiException("Leonardo is not healthy"));
    }
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo; this step only makes an api call to leonardo to check if its healthy
    return StepResult.getStepResultSuccess();
  }
}
