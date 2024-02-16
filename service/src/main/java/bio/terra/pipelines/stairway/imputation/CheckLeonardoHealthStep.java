package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;

/** This step queries the Leonardo status endpoint to check if it is healthy */
public class CheckLeonardoHealthStep implements Step {
  private final LeonardoService leonardoService;

  public CheckLeonardoHealthStep(LeonardoService leonardoService) {
    this.leonardoService = leonardoService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) throws RetryException {
    HealthCheck.Result healthResult = leonardoService.checkHealth();
    if (!healthResult.isOk()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_RETRY,
          new LeonardoServiceApiException(new ApiException("Leonardo is not healthy")));
    }
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) throws InterruptedException {
    // nothing to undo; this step only makes an api call to leonardo to check if its healthy
    return StepResult.getStepResultSuccess();
  }
}
