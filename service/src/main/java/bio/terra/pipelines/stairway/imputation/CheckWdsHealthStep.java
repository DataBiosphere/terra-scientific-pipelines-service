package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import org.databiosphere.workspacedata.client.ApiException;

/**
 * This step checks the health of the wds app associated with the passed in workspace id
 *
 * <p>this step expects wds uri to be provided in the working map
 */
public class CheckWdsHealthStep implements Step {
  private final WdsService wdsService;
  private final SamService samService;

  public CheckWdsHealthStep(WdsService wdsService, SamService samService) {
    this.wdsService = wdsService;
    this.samService = samService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {

    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationJobFlightMapKeys.WDS_URI);

    String wdsUri = workingMap.get(RunImputationJobFlightMapKeys.WDS_URI, String.class);

    HealthCheckWorkspaceApps.Result healthResult =
        wdsService.checkHealth(wdsUri, samService.getTspsServiceAccountToken());
    if (!healthResult.isOk()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_RETRY,
          new WdsServiceApiException(new ApiException("WDS is not healthy")));
    }
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) throws InterruptedException {
    // nothing to undo; this step just checks that wds is healthy
    return StepResult.getStepResultSuccess();
  }
}
