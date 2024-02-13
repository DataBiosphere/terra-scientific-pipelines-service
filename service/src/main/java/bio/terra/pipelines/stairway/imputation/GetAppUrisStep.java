package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;

/**
 * This step queries for the URI of the cbas and wds for a given workspace id. It then stores both
 * of the URIs in the working map for downstream steps.
 */
public class GetAppUrisStep implements Step {
  private final LeonardoService leonardoService;
  private final SamService samService;

  public GetAppUrisStep(LeonardoService leonardoService, SamService samService) {
    this.leonardoService = leonardoService;
    this.samService = samService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) throws RetryException {
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID);

    String controlWorkspaceId =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID, String.class);

    HealthCheck.Result healthResult = leonardoService.checkHealth();
    if (!healthResult.isOk()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_RETRY,
          new LeonardoServiceApiException(new ApiException("Leonardo is not healthy")));
    }

    List<ListAppResponse> appResponseList =
        leonardoService.getApps(controlWorkspaceId, samService.getTspsServiceAccountToken(), false);
    String cbasUri =
        leonardoService.getCbasUrlFromGetAppResponse(appResponseList, controlWorkspaceId);
    String wdsUri =
        leonardoService.getWdsUrlFromGetAppResponse(appResponseList, controlWorkspaceId);

    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(RunImputationJobFlightMapKeys.CBAS_URI, cbasUri);
    workingMap.put(RunImputationJobFlightMapKeys.WDS_URI, wdsUri);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) throws InterruptedException {
    // nothing to undo; this step only puts stuff in the working map for downstream steps
    return StepResult.getStepResultSuccess();
  }
}
