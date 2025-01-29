package bio.terra.pipelines.stairway.flights.imputation;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.steps.common.CompletePipelineRunStep;
import bio.terra.pipelines.stairway.steps.imputation.PrepareImputationInputsStep;
import bio.terra.pipelines.stairway.steps.imputation.azure.*;
import bio.terra.stairway.*;

public class RunImputationAzureJobFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 1, /* maxCount= */ 5);

  /** Retry for interacting with data plane apps */
  private final RetryRule dataPlaneAppRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 20, /* maxCount= */ 2);

  // addStep is protected in Flight, so make an override that is public
  @Override
  public void addStep(Step step, RetryRule retryRule) {
    super.addStep(step, retryRule);
  }

  public RunImputationAzureJobFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
    final FlightBeanBag flightBeanBag = FlightBeanBag.getFromObject(beanBag);

    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        JobMapKeys.PIPELINE_NAME,
        JobMapKeys.PIPELINE_ID,
        JobMapKeys.DOMAIN_NAME,
        JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK,
        JobMapKeys.DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK,
        JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK,
        ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
        ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
        ImputationJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        ImputationJobMapKeys.CONTROL_WORKSPACE_ID,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL,
        ImputationJobMapKeys.WDL_METHOD_NAME);

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(inputParameters.get(JobMapKeys.PIPELINE_NAME, String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    addStep(
        new PrepareImputationInputsStep(
            flightBeanBag.getPipelineInputsOutputsService(),
            flightBeanBag.getImputationConfiguration()),
        dbRetryRule);

    addStep(new CheckLeonardoHealthStep(flightBeanBag.getLeonardoService()), dataPlaneAppRetryRule);

    addStep(
        new GetAppUrisStep(flightBeanBag.getLeonardoService(), flightBeanBag.getSamService()),
        dataPlaneAppRetryRule);

    addStep(
        new CheckWdsHealthStep(flightBeanBag.getWdsService(), flightBeanBag.getSamService()),
        dataPlaneAppRetryRule);

    addStep(
        new AddWdsRowStep(flightBeanBag.getWdsService(), flightBeanBag.getSamService()),
        dataPlaneAppRetryRule);

    addStep(
        new CheckCbasHealthStep(flightBeanBag.getCbasService(), flightBeanBag.getSamService()),
        dataPlaneAppRetryRule);

    addStep(
        new SubmitCromwellRunSetStep(
            flightBeanBag.getCbasService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getPipelinesService(),
            flightBeanBag.getCbasConfiguration()),
        dataPlaneAppRetryRule);

    addStep(
        new PollCromwellRunSetStatusStep(
            flightBeanBag.getCbasService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getImputationConfiguration()),
        dataPlaneAppRetryRule);

    addStep(
        new FetchOutputsFromWdsStep(flightBeanBag.getWdsService(), flightBeanBag.getSamService()),
        dataPlaneAppRetryRule);

    addStep(new CompletePipelineRunStep(flightBeanBag.getPipelineRunsService()), dbRetryRule);
  }
}
