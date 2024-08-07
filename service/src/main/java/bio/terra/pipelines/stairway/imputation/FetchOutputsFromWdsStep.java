package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.databiosphere.workspacedata.model.RecordResponse;

/**
 * This step fetches outputs from a WDS record for a given job and stores them in the flight's
 * working map. These outputs are considered raw in that any cloud path values do not include SAS
 * tokens.
 */
public class FetchOutputsFromWdsStep implements Step {

  private final WdsService wdsService;
  private final SamService samService;

  public FetchOutputsFromWdsStep(WdsService wdsService, SamService samService) {
    this.wdsService = wdsService;
    this.samService = samService;
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
        RunImputationAzureJobFlightMapKeys.CONTROL_WORKSPACE_ID,
        RunImputationAzureJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS);

    String controlWorkspaceId =
        inputParameters.get(RunImputationAzureJobFlightMapKeys.CONTROL_WORKSPACE_ID, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);
    List<PipelineOutputDefinition> outputDefinitions =
        inputParameters.get(
            RunImputationAzureJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
            new TypeReference<>() {});

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationAzureJobFlightMapKeys.WDS_URI);

    String wdsUri = workingMap.get(RunImputationAzureJobFlightMapKeys.WDS_URI, String.class);

    RecordResponse recordResponse;
    try {
      recordResponse =
          wdsService.getRecord(
              wdsUri,
              samService.getTeaspoonsServiceAccountToken(),
              pipelineName.getValue(),
              controlWorkspaceId,
              jobId);
    } catch (WdsServiceException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    Map<String, String> outputs = new HashMap<>();
    for (PipelineOutputDefinition outputDefinition : outputDefinitions) {
      String keyName = outputDefinition.getName();
      String wdlVariableName = outputDefinition.getWdlVariableName();
      outputs.put(keyName, recordResponse.getAttributes().get(wdlVariableName).toString());
    }

    workingMap.put(RunImputationAzureJobFlightMapKeys.PIPELINE_RUN_OUTPUTS, outputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
