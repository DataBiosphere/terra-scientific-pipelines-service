package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.model.QuotaConsumptionRecord;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import io.sentry.Sentry;
import io.sentry.SentryLevel;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to mark a pipeline run as a success and write the outputs, effective quota consumed and
 * quota allocation deductions to the database
 *
 * <p>This step expects the JobMapKeys.USER_ID in the input parameters and
 * WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS and
 * WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS in the working map. The effective quota
 * consumed is the sum of the amounts in the quota consumption records.
 *
 * <p>WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE is an optional input in the working
 * map. If the output file sizes are present in the working map, it will write the file sizes to the
 * database, but if they are not present, it will continue without writing the file sizes since we
 * don't want to fail the entire step if we can't get the file sizes.
 */
public class CompletePipelineRunStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final QuotasService quotasService;
  private final Logger logger = LoggerFactory.getLogger(CompletePipelineRunStep.class);

  public CompletePipelineRunStep(
      PipelineRunsService pipelineRunsService, QuotasService quotasService) {
    this.pipelineRunsService = pipelineRunsService;
    this.quotasService = quotasService;
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
        WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS);
    Map<String, Object> outputsMap =
        workingMap.get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class);
    List<QuotaConsumptionRecord> quotaConsumptionRecords =
        workingMap.get(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, List.class);

    // calculate the total quota consumed for this run by summing the amounts in the quota
    // consumption records
    int quotaConsumed =
        quotaConsumptionRecords.stream().mapToInt(QuotaConsumptionRecord::amountConsumed).sum();

    // fetch output file sizes from working map, but if they are not present, continue
    // with an empty map since we don't want to fail the entire step if we can't get the
    // file sizes
    Map<String, Object> outputFileSizes = Collections.emptyMap();
    if (workingMap.containsKey(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE)) {
      outputFileSizes =
          workingMap.get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class);
    }

    PipelineRun pipelineRun =
        pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(
            jobId, userId, quotaConsumed, outputsMap, outputFileSizes);

    // record the quota deductions for this run in the database
    try {
      quotasService.recordQuotaDeductions(pipelineRun.getId(), quotaConsumptionRecords);
    } catch (RuntimeException e) {
      // the run and its quota consumption are already committed; don't fail a successful run
      // because writing deductions rows failed
      String msg =
          "Failed to record quota deductions for user %s, pipeline run %d and flight %s. Quota consumption records: %s"
              .formatted(
                  userId,
                  pipelineRun.getId(),
                  flightContext.getFlightId(),
                  quotaConsumptionRecords);
      logger.error(msg, e);
      Sentry.captureMessage(msg, SentryLevel.ERROR);
    }

    logger.info(
        "Marked run {} as a success and wrote job outputs and quota deductions to the db", jobId);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
