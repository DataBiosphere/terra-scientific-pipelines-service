package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.time.LocalDateTime;
import java.util.Map;

/**
 * This step creates or replaces a row to a workspace data table specific to the pipeline that was
 * launched. currently it writes the flight id as the primary key and all necessary wdl inputs
 * values.
 *
 * <p>this step expects the main tool config to be provided in the input parameter map
 *
 * <p>this step expects all pipeline inputs to be provided in the working map
 */
public class AddDataTableRowStep implements Step {
  private final RawlsService rawlsService;
  private final SamService samService;
  private final String toolConfigKey;

  public AddDataTableRowStep(
      RawlsService rawlsService, SamService samService, String toolConfigKey) {
    this.rawlsService = rawlsService;
    this.samService = samService;
    this.toolConfigKey = toolConfigKey;
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
        WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_NAME,
        toolConfigKey);

    String controlWorkspaceName =
        inputParameters.get(WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(
            WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    ToolConfig toolConfig = inputParameters.get(toolConfigKey, ToolConfig.class);
    String dataTableEntityName = toolConfig.dataTableEntityName();

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, WdlBasedPipelineJobMapKeys.ALL_PIPELINE_INPUTS);
    Map<String, Object> allPipelineInputs =
        workingMap.get(WdlBasedPipelineJobMapKeys.ALL_PIPELINE_INPUTS, new TypeReference<>() {});

    Entity entity =
        new Entity()
            .entityType(dataTableEntityName)
            .name(flightContext.getFlightId())
            .attributes(allPipelineInputs)
            .putAttributesItem("timestamp_start", LocalDateTime.now());
    try {
      rawlsService.upsertDataTableEntity(
          samService.getTeaspoonsServiceAccountToken(),
          controlWorkspaceProject,
          controlWorkspaceName,
          entity);
    } catch (RawlsServiceApiException e) {
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
