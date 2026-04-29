package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.Collections;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to mark a pipeline run as a success and write the outputs and effective quota consumed to
 * the database
 *
 * <p>This step expects the JobMapKeys.USER_ID in the input parameters and
 * WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS and
 * WdlBasedPipelineJobMapKeys.EFFECTIVE_QUOTA_CONSUMED in the working map.
 *
 * <p>WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE is an optional input in the working
 * map. If the output file sizes are present in the working map, it will write the file sizes to the
 * database, but if they are not present, it will continue without writing the file sizes since we
 * don't want to fail the entire step if we can't get the file sizes.
 */
public class CompletePipelineRunStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final Logger logger = LoggerFactory.getLogger(CompletePipelineRunStep.class);

  public CompletePipelineRunStep(PipelineRunsService pipelineRunsService) {
    this.pipelineRunsService = pipelineRunsService;
  }

  @Override
  @SuppressWarnings("java:S2259") // suppress warning for possible NPE when unboxing quotaConsumed,
  //  since we do validate that quotaConsumed is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID);

    UUID jobId = UUID.fromString(flightContext.getFlightId());
    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

    // validate and extract parameters from working map
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        workingMap,
        WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS,
        WdlBasedPipelineJobMapKeys.EFFECTIVE_QUOTA_CONSUMED);
    Map<String, String> outputsMap =
        workingMap.get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class);
    int quotaConsumed =
        workingMap.get(WdlBasedPipelineJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, Integer.class);

    // fetch output file sizes from working map, but if they are not present, continue
    // with an empty map since we don't want to fail the entire step if we can't get the
    // file sizes
    Map<String, Long> outputFileSizes = Collections.emptyMap();
    if (workingMap.containsKey(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE)) {
      outputFileSizes =
          workingMap.get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class);
    }

    pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(
        jobId, userId, quotaConsumed, outputsMap, outputFileSizes);

    logger.info("Marked run {} as a success and wrote job outputs & quota to the db", jobId);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
