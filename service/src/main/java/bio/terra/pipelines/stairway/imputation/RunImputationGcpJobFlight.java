package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.FetchQuotaConsumedFromDataTableStep;
import bio.terra.pipelines.stairway.PollQuotaConsumedSubmissionStatusStep;
import bio.terra.pipelines.stairway.SubmitQuotaConsumedSubmissionStep;
import bio.terra.pipelines.stairway.imputation.steps.CompletePipelineRunStep;
import bio.terra.pipelines.stairway.imputation.steps.PrepareImputationInputsStep;
import bio.terra.pipelines.stairway.imputation.steps.gcp.AddDataTableRowStep;
import bio.terra.pipelines.stairway.imputation.steps.gcp.FetchOutputsFromDataTableStep;
import bio.terra.pipelines.stairway.imputation.steps.gcp.PollCromwellSubmissionStatusStep;
import bio.terra.pipelines.stairway.imputation.steps.gcp.SubmitCromwellSubmissionStep;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.RetryRule;
import bio.terra.stairway.RetryRuleExponentialBackoff;
import bio.terra.stairway.RetryRuleFixedInterval;
import bio.terra.stairway.Step;

public class RunImputationGcpJobFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 1, /* maxCount= */ 5);

  /**
   * Use for a short exponential backoff retry, for operations that should be completable within a
   * few seconds.
   */
  private final RetryRule externalServiceRetryRule =
      // maxOperationTimeSeconds must be larger than socket timeout (20s), otherwise a socket
      // timeout
      // won't be retried.
      new RetryRuleExponentialBackoff(1, 8, /* maxOperationTimeSeconds */ 30);

  // addStep is protected in Flight, so make an override that is public
  @Override
  public void addStep(Step step, RetryRule retryRule) {
    super.addStep(step, retryRule);
  }

  public RunImputationGcpJobFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
    final FlightBeanBag flightBeanBag = FlightBeanBag.getFromObject(beanBag);

    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        JobMapKeys.PIPELINE_NAME,
        JobMapKeys.PIPELINE_ID,
        JobMapKeys.RESULT_PATH,
        JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK,
        JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK,
        ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
        ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
        ImputationJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL,
        ImputationJobMapKeys.WDL_METHOD_NAME,
        ImputationJobMapKeys.WDL_METHOD_VERSION);

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(inputParameters.get(JobMapKeys.PIPELINE_NAME, String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    addStep(
        new PrepareImputationInputsStep(
            flightBeanBag.getPipelinesService(), flightBeanBag.getImputationConfiguration()),
        dbRetryRule);

    addStep(
        new AddDataTableRowStep(flightBeanBag.getRawlsService(), flightBeanBag.getSamService()),
        externalServiceRetryRule);

    addStep(
        new SubmitQuotaConsumedSubmissionStep(
            flightBeanBag.getRawlsService(), flightBeanBag.getSamService()),
        externalServiceRetryRule);

    addStep(
        new PollQuotaConsumedSubmissionStatusStep(
            flightBeanBag.getRawlsService(), flightBeanBag.getSamService()),
        externalServiceRetryRule);

    addStep(
        new FetchQuotaConsumedFromDataTableStep(
            flightBeanBag.getRawlsService(), flightBeanBag.getSamService()),
        externalServiceRetryRule);

    addStep(
        new SubmitCromwellSubmissionStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getImputationConfiguration()),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getImputationConfiguration()),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService()),
        externalServiceRetryRule);

    addStep(new CompletePipelineRunStep(flightBeanBag.getPipelineRunsService()), dbRetryRule);
  }
}
