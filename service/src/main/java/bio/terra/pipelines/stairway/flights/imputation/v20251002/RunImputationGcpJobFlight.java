package bio.terra.pipelines.stairway.flights.imputation.v20251002;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.steps.common.*;
import bio.terra.pipelines.stairway.steps.imputation.AddDataTableRowStep;
import bio.terra.pipelines.stairway.steps.imputation.PrepareImputationInputsStep;
import bio.terra.stairway.*;

public class RunImputationGcpJobFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/* intervalSeconds= */ 1, /* maxCount= */ 5);

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
        JobMapKeys.PIPELINE_VERSION,
        JobMapKeys.PIPELINE_ID,
        JobMapKeys.DOMAIN_NAME,
        JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK,
        JobMapKeys.DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK,
        JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK,
        ImputationJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL,
        ImputationJobMapKeys.PIPELINE_TOOL_CONFIG,
        ImputationJobMapKeys.QUOTA_TOOL_CONFIG,
        ImputationJobMapKeys.INPUT_QC_TOOL_CONFIG);

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(inputParameters.get(JobMapKeys.PIPELINE_NAME, String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    Integer pipelineVersion = inputParameters.get(JobMapKeys.PIPELINE_VERSION, Integer.class);

    addStep(
        new PrepareImputationInputsStep(
            flightBeanBag.getPipelineInputsOutputsService(),
            flightBeanBag
                .getPipelineConfigurations()
                .getArrayImputation()
                .get(pipelineVersion.toString())),
        dbRetryRule);

    addStep(
        new AddDataTableRowStep(flightBeanBag.getRawlsService(), flightBeanBag.getSamService()),
        externalServiceRetryRule);

    // Check input for quota to be consumed
    addStep(
        new SubmitCromwellSubmissionStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            ImputationJobMapKeys.QUOTA_TOOL_CONFIG,
            ImputationJobMapKeys.QUOTA_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            ImputationJobMapKeys.QUOTA_TOOL_CONFIG,
            ImputationJobMapKeys.QUOTA_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService(),
            ImputationJobMapKeys.QUOTA_TOOL_CONFIG,
            ImputationJobMapKeys.QUOTA_OUTPUTS),
        externalServiceRetryRule);

    addStep(
        new QuotaConsumedValidationStep(
            flightBeanBag.getQuotasService(), flightBeanBag.getPipelineRunsService()),
        dbRetryRule);

    // run QC on user input
    addStep(
        new SubmitCromwellSubmissionStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            ImputationJobMapKeys.INPUT_QC_TOOL_CONFIG,
            ImputationJobMapKeys.INPUT_QC_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            ImputationJobMapKeys.INPUT_QC_TOOL_CONFIG,
            ImputationJobMapKeys.INPUT_QC_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService(),
            ImputationJobMapKeys.INPUT_QC_TOOL_CONFIG,
            ImputationJobMapKeys.INPUT_QC_OUTPUTS),
        externalServiceRetryRule);

    addStep(new InputQcValidationStep(), dbRetryRule);

    // run imputation
    addStep(
        new SubmitCromwellSubmissionStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            ImputationJobMapKeys.PIPELINE_TOOL_CONFIG,
            ImputationJobMapKeys.PIPELINE_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            ImputationJobMapKeys.PIPELINE_TOOL_CONFIG,
            ImputationJobMapKeys.PIPELINE_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService(),
            ImputationJobMapKeys.PIPELINE_TOOL_CONFIG,
            ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS),
        externalServiceRetryRule);

    addStep(new CompletePipelineRunStep(flightBeanBag.getPipelineRunsService()), dbRetryRule);

    addStep(
        new SendJobSucceededNotificationStep(flightBeanBag.getNotificationService()), dbRetryRule);
  }
}
