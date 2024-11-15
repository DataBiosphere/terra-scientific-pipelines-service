package bio.terra.pipelines.stairway;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.imputation.ImputationJobMapKeys;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step calls Rawls to fetch outputs from a data table row for a given quota consumed job. It
 * specifically fetches the quota consumed value from the data table row using the quota_consumed
 * key
 *
 * <p>This step expects nothing from the working map
 */
public class FetchQuotaConsumedFromDataTableStep implements Step {

  private final RawlsService rawlsService;
  private final SamService samService;
  private final Logger logger = LoggerFactory.getLogger(FetchQuotaConsumedFromDataTableStep.class);

  public FetchQuotaConsumedFromDataTableStep(RawlsService rawlsService, SamService samService) {
    this.rawlsService = rawlsService;
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
        JobMapKeys.PIPELINE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME);

    String controlWorkspaceBillingProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);

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

    // extract quota_consumed from entity
    int quotaConsumed;
    try {
      quotaConsumed = (int) entity.getAttributes().get("quota_consumed");
      if (quotaConsumed <= 0) {
        return new StepResult(
            StepStatus.STEP_RESULT_FAILURE_FATAL,
            new InternalServerErrorException("Quota consumed is unexpectedly not greater than 0"));
      }
    } catch (NullPointerException e) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException("Quota consumed is unexpectedly null"));
    }

    logger.info("Quota consumed: {}", quotaConsumed);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
