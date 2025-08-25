package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.*;
import java.util.List;
import java.util.Map;

/**
 * This step calls Rawls to fetch outputs from a data table row for a given quota consumed job. It
 * specifically fetches the quota consumed value from the data table row using the quota_consumed
 * key. If successful, it stores the quota consumed value in the working map.
 *
 * <p>This step expects nothing from the working map
 */
public class FetchValuesFromDataTableStep implements Step {

  private final RawlsService rawlsService;
  private final SamService samService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;

  private final String toolConfigKey;
  private final String toolOutputsKey;

  public FetchValuesFromDataTableStep(
      RawlsService rawlsService,
      SamService samService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      String toolConfigKey,
      String toolOutputsKey) {
    this.rawlsService = rawlsService;
    this.samService = samService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.toolConfigKey = toolConfigKey;
    this.toolOutputsKey = toolOutputsKey;
  }

  @Override
  @SuppressWarnings(
      "java:S2259") // suppress warning for possible NPE when calling pipelineName.getValue(),
  //  since we do validate that pipelineName is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) {
    String jobId = flightContext.getFlightId();

    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        toolConfigKey);

    String controlWorkspaceBillingProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
    ToolConfig toolConfig = inputParameters.get(toolConfigKey, ToolConfig.class);
    List<PipelineOutputDefinition> outputDefinitions = toolConfig.outputDefinitions();

    Entity entity;
    try {
      entity =
          rawlsService.getDataTableEntity(
              samService.getTeaspoonsServiceAccountToken(),
              controlWorkspaceBillingProject,
              controlWorkspaceName,
              pipelineName.getValue(),
              jobId);
    } catch (RawlsServiceApiException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    // this will throw an error and fail the task without retries if any of the output definitions
    // are missing or empty
    Map<String, String> outputs =
        pipelineInputsOutputsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);

    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(toolOutputsKey, outputs);

    //    // extract quota_consumed from entity
    //    int quotaConsumed;
    //    try {
    //      quotaConsumed = (int) entity.getAttributes().get("quota_consumed");
    //      if (quotaConsumed <= 0) {
    //        return new StepResult(
    //            StepStatus.STEP_RESULT_FAILURE_FATAL,
    //            new InternalServerErrorException("Quota consumed is unexpectedly not greater than
    // 0"));
    //      }
    //    } catch (NullPointerException e) {
    //      return new StepResult(
    //          StepStatus.STEP_RESULT_FAILURE_FATAL,
    //          new InternalServerErrorException("Quota consumed is unexpectedly null"));
    //    }
    //
    //    // store the raw quota consumed value in the working map to be used in a subsequent step
    //    workingMap.put(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, quotaConsumed);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
