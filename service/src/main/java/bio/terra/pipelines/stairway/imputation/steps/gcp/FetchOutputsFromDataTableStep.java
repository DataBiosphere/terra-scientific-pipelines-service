package bio.terra.pipelines.stairway.imputation.steps.gcp;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.List;
import java.util.Map;

/**
 * This step calls Rawls to fetch outputs from a data table row for a given job and stores them in
 * the flight's working map. These outputs are considered raw in that they are cloud paths and not
 * signed urls.
 */
public class FetchOutputsFromDataTableStep implements Step {

  private final RawlsService rawlsService;
  private final SamService samService;
  private final PipelineRunsService pipelineRunsService;

  public FetchOutputsFromDataTableStep(
      RawlsService rawlsService, SamService samService, PipelineRunsService pipelineRunsService) {
    this.rawlsService = rawlsService;
    this.samService = samService;
    this.pipelineRunsService = pipelineRunsService;
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
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME,
        RunImputationJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS);

    String controlWorkspaceBillingProject =
        inputParameters.get(
            RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    String controlWorkspaceName =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);
    List<PipelineOutputDefinition> outputDefinitions =
        inputParameters.get(
            RunImputationJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS, new TypeReference<>() {});

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

    // this will throw an error if any of the output definitions are not found or empty
    Map<String, String> outputs;
    try {
      outputs = pipelineRunsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);
    } catch (InternalServerErrorException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, e);
    }

    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(RunImputationJobFlightMapKeys.PIPELINE_RUN_OUTPUTS, outputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
