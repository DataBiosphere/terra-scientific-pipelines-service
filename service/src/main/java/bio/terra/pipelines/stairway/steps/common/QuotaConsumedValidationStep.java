package bio.terra.pipelines.stairway.steps.common;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.service.exception.PipelineCheckFailedException;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.*;
import java.util.Map;
import java.util.UUID;

/**
 * This step: 1. Extracts the raw quota detected from the QuotaConsumed wdl 2. Stores that
 * rawQuotaConsumed value in the PipelineRuns database for the pipeline run 3. Sets the
 * effectiveQuotaConsumed value to the maximum of the rawQuotaConsumed and the min_quota_consumed
 * for the pipeline being run 4. Checks that the effectiveQuotaConsumed for this run does not cause
 * the user to exceed their quota limit
 *
 * <p>If everything passes then this step writes the effective quota consumed for this run to the
 * working map and to the user quotas table.
 *
 * <p>This step expects the quota wdl outputs (map{quotaConsumed:value}) to be provided in the
 * working map
 */
public class QuotaConsumedValidationStep implements Step {
  private final QuotasService quotasService;
  private final PipelineRunsService pipelineRunsService;

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

    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters, JobMapKeys.PIPELINE_NAME, JobMapKeys.USER_ID);

    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.QUOTA_OUTPUTS);

    Map<?, ?> quotaOutputsMap = workingMap.get(ImputationJobMapKeys.QUOTA_OUTPUTS, Map.class);
    Object quotaConsumedObj = quotaOutputsMap.get("quotaConsumed");
    if (quotaConsumedObj == null) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException("Missing 'quotaConsumed' entry in quota outputs map."));
    }
    String rawQuotaConsumedValue = quotaConsumedObj.toString();
    int rawQuotaConsumed = Integer.parseInt(rawQuotaConsumedValue);
    if (rawQuotaConsumed <= 0) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException("Quota consumed is unexpectedly not greater than 0"));
    }

    // update the rawQuotaConsumed for this pipeline run in the db
    pipelineRunsService.setPipelineRunRawQuotaConsumed(
        UUID.fromString(flightContext.getFlightId()), userId, rawQuotaConsumed);

    // we want to have the ability to have a floor for quota consumed, so we take the max of the
    // pipeline's min quota consumed and the quota consumed for this run
    Integer quotaUsedForThisRun =
        Math.max(
            quotasService.getPipelineQuota(pipelineName).getMinQuotaConsumed(), rawQuotaConsumed);

    // check if user quota used plus quota consumed is less than or equal to user quota
    UserQuota userQuota = quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipelineName);

    // user quota has been exceeded, fail the flight
    int totalQuotaConsumed = userQuota.getQuotaConsumed() + quotaUsedForThisRun;
    if (totalQuotaConsumed > userQuota.getQuota()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new PipelineCheckFailedException(
              String.format(
                  "User quota exceeded for pipeline %s. User quota limit: %d, Quota consumed before this run: %d, "
                      + "Quota consumed for this run: %d.  If you would like to request a quota increase, you can "
                      + "email scientific-services-support@broadinstitute.org",
                  pipelineName.getValue(),
                  userQuota.getQuota(),
                  userQuota.getQuotaConsumed(),
                  quotaUsedForThisRun)));
    }
    // quota has not been exceeded, update user quota consumed
    quotasService.updateQuotaConsumed(userQuota, totalQuotaConsumed);

    // store the raw and effective quota consumed values in the working map
    workingMap.put(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, rawQuotaConsumed);
    workingMap.put(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, quotaUsedForThisRun);

    return StepResult.getStepResultSuccess();
  }

  // Undo the increase in quota consumed
  @SuppressWarnings("java:S2259") // suppress warning for possible NPE when calling workingMap.get,
  // these values should exist if this step succeeded and is getting undone due to a
  // future step failing.
  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // grab variable needed
    var inputParameters = flightContext.getInputParameters();
    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    FlightMap workingMap = flightContext.getWorkingMap();

    // update the user quota to be what it was before this run's quota was added
    Integer quotaUsedForThisRun =
        workingMap.get(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, Integer.class);

    UserQuota userQuota = quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipelineName);

    quotasService.updateQuotaConsumed(
        userQuota, userQuota.getQuotaConsumed() - quotaUsedForThisRun);
    return StepResult.getStepResultSuccess();
  }
}
