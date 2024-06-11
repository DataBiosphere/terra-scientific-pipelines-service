package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
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
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID);

    String controlWorkspaceId =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationJobFlightMapKeys.WDS_URI);

    String wdsUri = workingMap.get(RunImputationJobFlightMapKeys.WDS_URI, String.class);

    RecordResponse recordResponse;
    try {
      recordResponse =
          wdsService.getRecord(
              wdsUri,
              samService.getTspsServiceAccountToken(),
              pipelineName.getValue(),
              controlWorkspaceId,
              jobId);
    } catch (WdsServiceException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    // in TSPS-197 we will get the outputKeys from the pipelines db table
    List<String> outputKeys =
        List.of("imputed_multi_sample_vcf", "imputed_multi_sample_vcf_index", "chunks_info");

    Map<String, String> outputs =
        recordResponse.getAttributes().entrySet().stream()
            .filter(entry -> outputKeys.contains(entry.getKey()))
            .collect(Collectors.toMap(Map.Entry::getKey, entry -> entry.getValue().toString()));

    workingMap.put(RunImputationJobFlightMapKeys.RAW_OUTPUTS_MAP, outputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
