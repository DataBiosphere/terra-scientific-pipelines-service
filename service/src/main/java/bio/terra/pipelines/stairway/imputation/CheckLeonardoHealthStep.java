package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;

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
          new LeonardoServiceApiException("Leonardo is not healthy"));
    }
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // this is the first step in RunImputationJobFlight.
    // increment pipeline failed counter if undoStep is called which means the flight failed
    // to be moved to a StairwayHook in https://broadworkbench.atlassian.net/browse/TSPS-181
    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(
            flightContext
                .getInputParameters()
                .get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    MetricsUtils.incrementPipelineRunFailed(pipelinesEnum);

    // nothing to undo; this step only makes an api call to leonardo to check if its healthy
    return StepResult.getStepResultSuccess();
  }
}
