package bio.terra.pipelines.stairway.steps.imputation.azure;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.*;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;

/**
 * This step queries for the URI of the cbas and wds for a given workspace id. It then stores both
 * of the URIs in the working map for downstream steps.
 *
 * <p>this step expects control workspace id to provided in the input parameter map
 *
 * <p>this step writes cbas uri and wds uri to the working map
 */
public class GetAppUrisStep implements Step {
  private final LeonardoService leonardoService;
  private final SamService samService;

  public GetAppUrisStep(LeonardoService leonardoService, SamService samService) {
    this.leonardoService = leonardoService;
    this.samService = samService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, ImputationJobMapKeys.CONTROL_WORKSPACE_ID);

    String controlWorkspaceId =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_ID, String.class);

    List<ListAppResponse> appResponseList =
        leonardoService.getApps(controlWorkspaceId, samService.getTeaspoonsServiceAccountToken());
    String cbasUri =
        leonardoService.getCbasUrlFromGetAppResponse(appResponseList, controlWorkspaceId);
    String wdsUri =
        leonardoService.getWdsUrlFromGetAppResponse(appResponseList, controlWorkspaceId);

    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(ImputationJobMapKeys.CBAS_URI, cbasUri);
    workingMap.put(ImputationJobMapKeys.WDS_URI, wdsUri);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; this step only puts stuff in the working map for downstream steps
    return StepResult.getStepResultSuccess();
  }
}
