package bio.terra.pipelines.stairway;

import bio.terra.common.exception.BadRequestException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.imputation.ImputationJobMapKeys;
import bio.terra.stairway.*;

/**
 * This step calls Rawls to fetch outputs from a data table row for a given quota consumed job. It
 * specifically fetches the quota consumed value from the data table row using the quota_consumed
 * key. If successful, it stores the quota consumed value in the working map.
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
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.QUOTA_CONSUMED);

    // check if user quota used plus quota consumed is less than or equal to user quota
    Integer quotaUsedForThisRun =
        workingMap.get(ImputationJobMapKeys.QUOTA_CONSUMED, Integer.class);
    UserQuota userQuota = quotasService.getQuotaForUserAndPipeline(userId, pipelineName);

    // user quota has been exceeded, fail the flight
    int totalQuotaConsumed = userQuota.getQuotaConsumed() + quotaUsedForThisRun;
    if (totalQuotaConsumed > userQuota.getQuota()) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new BadRequestException("User quota exceeded for pipeline " + pipelineName.getValue()));
    }
    // quota has not been exceeded, update user quota consumed
    quotasService.updateQuotaConsumed(userQuota, totalQuotaConsumed);

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
    int quotaUsedForThisRun = workingMap.get(ImputationJobMapKeys.QUOTA_CONSUMED, Integer.class);

    UserQuota userQuota = quotasService.getQuotaForUserAndPipeline(userId, pipelineName);

    quotasService.updateQuotaConsumed(
        userQuota, userQuota.getQuotaConsumed() - quotaUsedForThisRun);
    return StepResult.getStepResultSuccess();
  }
}
