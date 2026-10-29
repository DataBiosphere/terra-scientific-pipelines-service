package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.enumFromPipelineKey;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.model.PipelineQuota;
import bio.terra.pipelines.model.QuotaConsumptionRecord;
import bio.terra.pipelines.model.UserQuotaTotals;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.service.exception.PipelineCheckFailedException;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.pipelines.stairway.steps.exception.InternalStepException;
import bio.terra.stairway.*;
import io.sentry.Sentry;
import io.sentry.SentryLevel;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.dao.TransientDataAccessException;

/**
 * This step:
 *
 * <p>1. Extracts the raw quota detected from the QuotaConsumed wdl
 *
 * <p>2. Stores that rawQuotaConsumed value in the PipelineRuns database for the pipeline run
 *
 * <p>3. Sets the effectiveQuotaConsumed value to the maximum of the rawQuotaConsumed and the
 * min_quota_consumed for the pipeline being run
 *
 * <p>4. Checks that the effectiveQuotaConsumed for this run does not cause the user to exceed their
 * quota limit
 *
 * <p>If everything passes then this step consumes the effectiveQuotaConsumed from the user's quota
 * allocations for this pipeline and stores the QuotaConsumptionRecords in the working map.
 *
 * <p>This step expects the quota wdl outputs (map{quotaConsumed:value}) to be provided in the
 * working map
 */
public class QuotaConsumedValidationStep implements Step {
  private final QuotasService quotasService;
  private final PipelineRunsService pipelineRunsService;

  private static final Logger logger = LoggerFactory.getLogger(QuotaConsumedValidationStep.class);

  public QuotaConsumedValidationStep(
      QuotasService quotasService, PipelineRunsService pipelineRunsService) {
    this.quotasService = quotasService;
    this.pipelineRunsService = pipelineRunsService;
  }

  @Override
  @SuppressWarnings(
      "java:S2259") // suppress warning for possible NPE when calling pipelineName.getValue(),
  //  since we do validate that pipelineName is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) {
    String flightId = flightContext.getFlightId();

    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters, JobMapKeys.PIPELINE_KEY, JobMapKeys.USER_ID);

    PipelinesEnum pipelineName =
        enumFromPipelineKey(inputParameters.get(JobMapKeys.PIPELINE_KEY, String.class));
    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS);

    Map<?, ?> quotaOutputsMap = workingMap.get(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, Map.class);
    Object quotaConsumedObj = quotaOutputsMap.get("quotaConsumed");
    if (quotaConsumedObj == null) {
      logger.error("Missing 'quotaConsumed' entry in quota outputs map for flight {}.", flightId);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalStepException("Quota consumption information not found."));
    }
    String rawQuotaConsumedValue = quotaConsumedObj.toString();
    int rawQuotaConsumed = Integer.parseInt(rawQuotaConsumedValue);
    if (rawQuotaConsumed <= 0) {
      logger.error(
          "Quota consumed is unexpectedly not greater than 0 for flight {}. Quota consumed value: {}",
          flightId,
          rawQuotaConsumed);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalStepException("No valid input files found; please check your input."));
    }

    PipelineQuota pipelineQuota = quotasService.getPipelineQuota(pipelineName);

    if (rawQuotaConsumed > pipelineQuota.getMaxQuotaConsumed()) {
      logger.error(
          "Quota consumed value of {} exceeds the maximum quota consumed of {} for pipeline {} for flight {}.",
          rawQuotaConsumed,
          pipelineQuota.getMaxQuotaConsumed(),
          pipelineName.getConfigKeyValue(),
          flightId);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new PipelineCheckFailedException(
              String.format(
                  "The number of submitted %s for this job (%d) exceeds the maximum allowed (%d). Please try again with a smaller input size.",
                  pipelineQuota.getQuotaUnits().getValue(),
                  rawQuotaConsumed,
                  pipelineQuota.getMaxQuotaConsumed())));
    }

    // update the rawQuotaConsumed for this pipeline run in the db
    pipelineRunsService.setPipelineRunRawQuotaConsumed(
        UUID.fromString(flightContext.getFlightId()), userId, rawQuotaConsumed);

    // we want to have the ability to have a floor for quota consumed, so we take the max of the
    // pipeline's min quota consumed and the quota consumed for this run
    Integer quotaUsedForThisRun = Math.max(pipelineQuota.getMinQuotaConsumed(), rawQuotaConsumed);

    // check if the user has enough quota remaining for this run
    UserQuotaTotals userQuotaTotals =
        quotasService.getQuotaTotalsForUserAndPipeline(userId, pipelineName);
    int totalQuotaConsumed = userQuotaTotals.totalConsumed() + quotaUsedForThisRun;

    // helper message for quota exceeded error, to be used in both the totals check and the consume
    // operation
    String quotaExceededMsg =
        String.format(
            "User quota exceeded for pipeline %s. User quota limit: %d, Quota consumed before this run: %d, "
                + "Quota consumed for this run: %d. Please visit https://services.terra.bio/pipelines/quotas to purchase additional quota.",
            pipelineName.getLowerCaseValue(),
            userQuotaTotals.totalAllocated(),
            userQuotaTotals.totalConsumed(),
            quotaUsedForThisRun);

    // if the total quota consumed exceeds the user's allocated quota, fail the step
    if (totalQuotaConsumed > userQuotaTotals.totalAllocated()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new PipelineCheckFailedException(quotaExceededMsg));
    }

    // update the quota for this run from the user's allocations for this pipeline
    List<QuotaConsumptionRecord> quotaConsumptionRecords;
    try {
      quotaConsumptionRecords =
          quotasService.consumeQuotaFromAllocationsForUserAndPipeline(
              userId, pipelineName, quotaUsedForThisRun);
    } catch (BadRequestException e) {
      // quota ran out between the totals check and the consume operation (e.g. concurrent job)
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new PipelineCheckFailedException(quotaExceededMsg));
    } catch (TransientDataAccessException e) {
      logger.warn("Transient DB error consuming quota for flight {}. Retrying.", flightId, e);
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    } catch (InternalServerErrorException e) {
      // consumeQuotaFromAllocationsForUserAndPipeline can throw this if there is a data integrity
      // issue in an allocation record, which should never happen, but if it does, we want to fail
      // the step and not proceed with the pipeline run
      logger.error("Failed to consume quota for flight {}", flightId, e);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalStepException("Error updating quota for job."));
    }

    // if no quota consumption records were returned, fail the step
    // this should never happen, but if it does, we want to fail the step and not proceed with the
    // pipeline run
    if (quotaConsumptionRecords == null || quotaConsumptionRecords.isEmpty()) {
      logger.error(
          "No quota consumption records returned for user {} and pipeline {} for flight {}.",
          userId,
          pipelineName.getLowerCaseValue(),
          flightId);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalStepException(
              "Something went wrong while updating quota for user %s and pipeline %s for flight %s."
                  .formatted(userId, pipelineName.getLowerCaseValue(), flightId)));
    }

    // store the raw quota consumed and the quota consumption records in the working map
    workingMap.put(WdlBasedPipelineJobMapKeys.RAW_QUOTA_CONSUMED, rawQuotaConsumed);
    workingMap.put(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, quotaConsumptionRecords);

    return StepResult.getStepResultSuccess();
  }

  // Undo the quota consumed for this job run
  @SuppressWarnings({
    "java:S2259", // suppress warning for possible NPE when calling workingMap.get,
    // these values should exist if this step succeeded and is getting undone due to a
    // future step failing.
    "unchecked" // the records are stored as a List<QuotaConsumptionRecord>, but read back as a raw
    // List because Stairway's mapper writes type ids for each element
  })
  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // grab variable needed
    var inputParameters = flightContext.getInputParameters();
    PipelinesEnum pipelineName =
        enumFromPipelineKey(inputParameters.get(JobMapKeys.PIPELINE_KEY, String.class));

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    FlightMap workingMap = flightContext.getWorkingMap();

    List<QuotaConsumptionRecord> records =
        workingMap.get(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, List.class);

    if (records == null || records.isEmpty()) {
      logger.info(
          "No quota consumption records found for user {}, pipeline {} and flight {}. Nothing to undo.",
          userId,
          pipelineName.getLowerCaseValue(),
          flightContext.getFlightId());
      return StepResult.getStepResultSuccess();
    }

    // revert the quota consumed for this job run
    try {
      quotasService.revertQuotaConsumption(records);
    } catch (TransientDataAccessException e) {
      // the revert runs in a single transaction that has been rolled back, so it is safe to retry
      logger.warn(
          "Transient DB error reverting quota for flight {}. Retrying.",
          flightContext.getFlightId(),
          e);
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    } catch (RuntimeException e) {
      // revertQuotaConsumption can throw if an allocation is missing, if reverting would make
      // consumed quota negative, or for unexpected DB errors. In all cases the quota consumed by
      // this job has not been returned, so log the records needed to fix it manually and alert.
      String errorMessage =
          "Failed to revert quota consumption for user %s, pipeline %s and flight %s. Quota consumption records: %s"
              .formatted(
                  userId, pipelineName.getLowerCaseValue(), flightContext.getFlightId(), records);
      logger.error(errorMessage, e);
      Sentry.captureMessage(errorMessage, SentryLevel.ERROR);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalStepException("Error while updating quota for job."));
    }

    logger.info(
        "Reverted quota consumption for user {}, pipeline {} and flight {}.",
        userId,
        pipelineName.getLowerCaseValue(),
        flightContext.getFlightId());

    return StepResult.getStepResultSuccess();
  }
}
