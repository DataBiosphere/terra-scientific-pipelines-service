package bio.terra.pipelines.stairway.imputation.gcp;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.*;

/**
 * This step creates or replaces a row to a workspace data table specific to the pipeline that was
 * launched. currently it writes the flight id as the primary key and a hardcoded "scatter" value
 * that will be replaced once inputs are being passed in from the user.
 *
 * <p>this step expects pipeline name and control workspace id to provided in the input parameter
 * map
 *
 * <p>this step expects wds uri to be provided in the working map
 */
public class AddDataTableRowStep implements Step {
  private final RawlsService rawlsService;
  private final SamService samService;

  public AddDataTableRowStep(RawlsService rawlsService, SamService samService) {
    this.rawlsService = rawlsService;
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
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_PROJECT,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME);

    String controlWorkspaceName =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_PROJECT, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    Entity entity =
        new Entity()
            .entityType(pipelineName.getValue())
            .name(flightContext.getFlightId())
            .putAttributesItem(
                "multisample_arrays_vcf",
                "gs://fc-secure-972057b6-fce1-4ec9-be0c-0e5c97fdc3e8/input_data/merged_500_input_samples.liftedover.vcf.gz");
    rawlsService.upsertDataTableEntity(
        samService.getTeaspoonsServiceAccountToken(),
        controlWorkspaceProject,
        controlWorkspaceName,
        entity);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; we don't need to remove the row that was added to WDS as it could be useful
    // for debugging. this may change in the future
    return StepResult.getStepResultSuccess();
  }
}
