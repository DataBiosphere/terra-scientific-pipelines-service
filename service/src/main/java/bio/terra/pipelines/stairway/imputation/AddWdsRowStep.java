package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import org.databiosphere.workspacedata.client.ApiException;
import org.databiosphere.workspacedata.model.RecordAttributes;
import org.databiosphere.workspacedata.model.RecordRequest;

/**
 * This step creates or replaces a row to a WDS table specific to the pipeline that was launched
 * currently it writes the flight id as the primary key and a hardcoded "scatter" value that will be
 * replaced once inputs are being passed in from the user.
 *
 * <p>The step expects a cbas uri to be passed in through the working map
 */
public class AddWdsRowStep implements Step {
  private final WdsService wdsService;
  private final SamService samService;

  public AddWdsRowStep(WdsService wdsService, SamService samService) {
    this.wdsService = wdsService;
    this.samService = samService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    FlightMap inputParameters = flightContext.getInputParameters();

    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID);

    String controlWorkspaceId =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID, String.class);
    String pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class);

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

    // hardcoded for now until we are using inputs from user
    RecordAttributes recordAttributes = new RecordAttributes();
    recordAttributes.put("scatter", 2);
    RecordRequest createRecordRequest = new RecordRequest().attributes(recordAttributes);
    try {
      wdsService.createOrReplaceRecord(
          wdsUri,
          samService.getTspsServiceAccountToken(),
          createRecordRequest,
          controlWorkspaceId,
          pipelineName,
          flightContext.getFlightId(),
          "flight_id");
    } catch (WdsServiceException e) {
      // not sure what exception makes sense to throw here
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, e);
    }

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) throws InterruptedException {
    // nothing to undo; we don't need to remove the row that was added to WDS as it could be useful
    // for debugging. this may change in the future
    return StepResult.getStepResultSuccess();
  }
}
