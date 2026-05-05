package bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.v20260428;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.pipelines.stairway.steps.common.AddDataTableRowStep;
import bio.terra.pipelines.stairway.steps.common.CompletePipelineRunStep;
import bio.terra.pipelines.stairway.steps.common.FetchOutputsFromDataTableStep;
import bio.terra.pipelines.stairway.steps.common.InputQcValidationStep;
import bio.terra.pipelines.stairway.steps.common.PollCromwellSubmissionStatusStep;
import bio.terra.pipelines.stairway.steps.common.PopulateFileOutputSizeStep;
import bio.terra.pipelines.stairway.steps.common.PrepareInputsStep;
import bio.terra.pipelines.stairway.steps.common.QuotaConsumedValidationStep;
import bio.terra.pipelines.stairway.steps.common.SendJobSucceededNotificationStep;
import bio.terra.pipelines.stairway.steps.common.SubmitCromwellSubmissionStep;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.RetryRule;
import bio.terra.stairway.RetryRuleExponentialBackoff;
import bio.terra.stairway.RetryRuleFixedInterval;
import bio.terra.stairway.Step;

public class RunWdlBasedPipelineJobFlight extends Flight {

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

  public RunWdlBasedPipelineJobFlight(FlightMap inputParameters, Object beanBag) {
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
        WdlBasedPipelineJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_NAME,
        WdlBasedPipelineJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        WdlBasedPipelineJobMapKeys.PIPELINE_TOOL_CONFIG,
        WdlBasedPipelineJobMapKeys.QUOTA_TOOL_CONFIG,
        WdlBasedPipelineJobMapKeys.INPUT_QC_TOOL_CONFIG);

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(inputParameters.get(JobMapKeys.PIPELINE_NAME, String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    Integer pipelineVersion = inputParameters.get(JobMapKeys.PIPELINE_VERSION, Integer.class);

    // prepare inputs is custom to pipeline
    PipelineConfigurations.WdlBasedPipelineConfig wdlBasedPipelineConfig;
    if (pipelinesEnum.equals(PipelinesEnum.ARRAY_IMPUTATION)) {
      wdlBasedPipelineConfig =
          flightBeanBag
              .getPipelineConfigurations()
              .getArrayImputation()
              .get(pipelineVersion.toString());
    } else {
      throw new IllegalArgumentException(
          String.format("Unsupported pipeline %s", pipelinesEnum.name()));
    }

    addStep(
        new PrepareInputsStep(
            flightBeanBag.getPipelineInputsOutputsService(), wdlBasedPipelineConfig),
        dbRetryRule);

    addStep(
        new AddDataTableRowStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            WdlBasedPipelineJobMapKeys.PIPELINE_TOOL_CONFIG),
        externalServiceRetryRule);

    // Check input for quota to be consumed
    addStep(
        new SubmitCromwellSubmissionStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            WdlBasedPipelineJobMapKeys.QUOTA_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.QUOTA_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            WdlBasedPipelineJobMapKeys.QUOTA_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.QUOTA_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService(),
            WdlBasedPipelineJobMapKeys.QUOTA_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS),
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
            WdlBasedPipelineJobMapKeys.INPUT_QC_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.INPUT_QC_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            WdlBasedPipelineJobMapKeys.INPUT_QC_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.INPUT_QC_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService(),
            WdlBasedPipelineJobMapKeys.INPUT_QC_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.INPUT_QC_OUTPUTS),
        externalServiceRetryRule);

    addStep(new InputQcValidationStep(), dbRetryRule);

    // run imputation
    addStep(
        new SubmitCromwellSubmissionStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            WdlBasedPipelineJobMapKeys.PIPELINE_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.PIPELINE_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            WdlBasedPipelineJobMapKeys.PIPELINE_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.PIPELINE_SUBMISSION_ID),
        externalServiceRetryRule);

    addStep(
        new FetchOutputsFromDataTableStep(
            flightBeanBag.getRawlsService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelineInputsOutputsService(),
            WdlBasedPipelineJobMapKeys.PIPELINE_TOOL_CONFIG,
            WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS),
        externalServiceRetryRule);

    // populate file sizes for pipeline outputs
    addStep(
        new PopulateFileOutputSizeStep(
            flightBeanBag.getPipelinesService(), flightBeanBag.getPipelineInputsOutputsService()),
        externalServiceRetryRule);

    addStep(new CompletePipelineRunStep(flightBeanBag.getPipelineRunsService()), dbRetryRule);

    addStep(
        new SendJobSucceededNotificationStep(flightBeanBag.getNotificationService()), dbRetryRule);
  }
}
