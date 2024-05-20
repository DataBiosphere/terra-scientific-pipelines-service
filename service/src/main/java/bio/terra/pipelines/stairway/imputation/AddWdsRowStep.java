package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
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
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID);

    String controlWorkspaceId =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationJobFlightMapKeys.WDS_URI);

    String wdsUri = workingMap.get(RunImputationJobFlightMapKeys.WDS_URI, String.class);

    // hardcoded for now until we are using inputs from user.
    // we are using the flight id as the primary key in the table created in WDS
    RecordAttributes recordAttributes = new RecordAttributes();
    recordAttributes.put("multi_sample_vcf", "a_fake_file.vcf.gz");
    recordAttributes.put("output_basename", "palantir_42_samples.hg38");
    RecordRequest createRecordRequest = new RecordRequest().attributes(recordAttributes);
    try {
      wdsService.createOrReplaceRecord(
          wdsUri,
          samService.getTspsServiceAccountToken(),
          createRecordRequest,
          controlWorkspaceId,
          pipelineName.getValue(),
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
