package bio.terra.pipelines.stairway.steps.common;

import bio.terra.common.exception.BadRequestException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.*;

/**
 * This step checks that the quota consumed for the flight is at least the mim_quota_consumed for
 * the pipeline being run. Once that is evaluated it then chekcs that the quota consumed for this
 * run does not cause the user to exceed their quota limit. If everything passes then this step
 * writes the effective quota consumed for this run to the working map.
 *
 * <p>This step expects quota consumed to be provided in the working map
 */
public class QuotaConsumedValidationStep implements Step {
  private final QuotasService quotasService;

  public QuotaConsumedValidationStep(QuotasService quotasService) {
    this.quotasService = quotasService;
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
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.RAW_QUOTA_CONSUMED);

    // we want to have the ability to have a floor for quota consumed, so we take the max of the
    // pipeline's min quota consumed and the quota consumed for this run
    Integer quotaUsedForThisRun =
        Math.max(
            quotasService.getPipelineQuota(pipelineName).getMinQuotaConsumed(),
            workingMap.get(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, Integer.class));

    // check if user quota used plus quota consumed is less than or equal to user quota
    UserQuota userQuota = quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipelineName);

    // user quota has been exceeded, fail the flight
    int totalQuotaConsumed = userQuota.getQuotaConsumed() + quotaUsedForThisRun;
    if (totalQuotaConsumed > userQuota.getQuota()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new BadRequestException(
              String.format(
                  "User quota exceeded for pipeline %s. User quota limit: %d, Quota consumed before this run: %d, "
                      + "Quota consumed for this run: %d.  If you would like to request a quota increase, you can "
                      + "email teaspoons-developers@broadinstitute.org",
                  pipelineName.getValue(),
                  userQuota.getQuota(),
                  userQuota.getQuotaConsumed(),
                  quotaUsedForThisRun)));
    }
    // quota has not been exceeded, update user quota consumed
    quotasService.updateQuotaConsumed(userQuota, totalQuotaConsumed);

    // store the effective quota consumed value in the working map to used in a subsequent step
    workingMap.put(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, quotaUsedForThisRun);

    return StepResult.getStepResultSuccess();
  }

  // Undo the increase in quota consumed
  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // grab variable needed
    var inputParameters = flightContext.getInputParameters();
    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    FlightMap workingMap = flightContext.getWorkingMap();

    // update the user quota to be what it was before this run's quota was added
    int quotaUsedForThisRun =
        workingMap.get(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, Integer.class);

    UserQuota userQuota = quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipelineName);

    quotasService.updateQuotaConsumed(
        userQuota, userQuota.getQuotaConsumed() - quotaUsedForThisRun);
    return StepResult.getStepResultSuccess();
  }
}
