package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Map;
import org.databiosphere.workspacedata.model.RecordAttributes;
import org.databiosphere.workspacedata.model.RecordRequest;

/**
 * This step creates or replaces a row to a WDS table specific to the pipeline that was launched
 * currently it writes the flight id as the primary key and a hardcoded "scatter" value that will be
 * replaced once inputs are being passed in from the user.
 *
 * <p>this step expects pipeline name and control workspace id to provided in the input parameter
 * map
 *
 * <p>this step expects wds uri to be provided in the working map
 */
public class AddWdsRowStep implements Step {
  private final WdsService wdsService;
  private final SamService samService;

  public AddWdsRowStep(WdsService wdsService, SamService samService) {
    this.wdsService = wdsService;
    this.samService = samService;
  }

  @Override
  @SuppressWarnings(
      "java:S2259") // suppress warning for possible NPE when calling pipelineName.getValue(),
  //  since we do validate that pipelineName is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationAzureJobFlightMapKeys.CONTROL_WORKSPACE_ID);

    String controlWorkspaceId =
        inputParameters.get(RunImputationAzureJobFlightMapKeys.CONTROL_WORKSPACE_ID, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        workingMap,
        RunImputationAzureJobFlightMapKeys.WDS_URI,
        RunImputationAzureJobFlightMapKeys.ALL_PIPELINE_INPUTS);

    String wdsUri = workingMap.get(RunImputationAzureJobFlightMapKeys.WDS_URI, String.class);
    Map<String, Object> allPipelineInputs =
        workingMap.get(
            RunImputationAzureJobFlightMapKeys.ALL_PIPELINE_INPUTS, new TypeReference<>() {});

    // create row to write to WDS
    RecordAttributes recordAttributes = new RecordAttributes();
    recordAttributes.putAll(allPipelineInputs);

    // add a timestamp - TSPS-227 will include generating a real timestamp that we can use here
    // instead
    recordAttributes.put("timestamp_start", System.currentTimeMillis());

    RecordRequest createRecordRequest = new RecordRequest().attributes(recordAttributes);
    try {
      wdsService.createOrReplaceRecord(
          wdsUri,
          samService.getTeaspoonsServiceAccountToken(),
          createRecordRequest,
          controlWorkspaceId,
          pipelineName.getValue(),
          flightContext.getFlightId(), // this is the primary key for WDS
          "flight_id");
    } catch (WdsServiceException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; we don't need to remove the row that was added to WDS as it could be useful
    // for debugging. this may change in the future
    return StepResult.getStepResultSuccess();
  }
}
