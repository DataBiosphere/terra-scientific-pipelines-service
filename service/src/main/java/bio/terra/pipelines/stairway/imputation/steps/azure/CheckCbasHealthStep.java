package bio.terra.pipelines.stairway.imputation.steps.azure;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.stairway.*;

/**
 * This step checks the health of the cbas app associated with the passed in workspace id
 *
 * <p>This step expects the cbas uri to be provided in the working map
 */
public class CheckCbasHealthStep implements Step {
  private final CbasService cbasService;
  private final SamService samService;

  public CheckCbasHealthStep(CbasService cbasService, SamService samService) {
    this.cbasService = cbasService;
    this.samService = samService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationJobFlightMapKeys.CBAS_URI);

    String cbasUri = workingMap.get(RunImputationJobFlightMapKeys.CBAS_URI, String.class);

    HealthCheckWorkspaceApps.Result healthResult =
        cbasService.checkHealth(cbasUri, samService.getTeaspoonsServiceAccountToken());

    if (!healthResult.isOk()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_RETRY, new CbasServiceApiException("CBAS is not healthy"));
    }

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; this step just queries for the health of a cbas app
    return StepResult.getStepResultSuccess();
  }
}
