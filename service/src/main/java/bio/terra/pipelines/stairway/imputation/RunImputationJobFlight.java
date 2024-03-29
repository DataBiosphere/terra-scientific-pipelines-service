package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.*;
import bio.terra.stairway.*;

public class RunImputationJobFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 1, /* maxCount= */ 5);

  /** Retry for interacting with data plane apps */
  private final RetryRule dataPlaneAppRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 20, /* maxCount= */ 5);

  // addStep is protected in Flight, so make an override that is public
  @Override
  public void addStep(Step step, RetryRule retryRule) {
    super.addStep(step, retryRule);
  }

  public RunImputationJobFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
    final FlightBeanBag flightBeanBag = FlightBeanBag.getFromObject(beanBag);

    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID.getKeyName(),
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_ID,
        RunImputationJobFlightMapKeys.PIPELINE_INPUTS,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID,
        RunImputationJobFlightMapKeys.WDL_METHOD_NAME,
        JobMapKeys.RESULT_PATH.getKeyName());

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    // write the job metadata to the Jobs table
    addStep(new WriteJobToDbStep(flightBeanBag.getImputationService()), dbRetryRule);

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
        new SubmitCromwellRunSetStep(flightBeanBag.getCbasService(), flightBeanBag.getSamService()),
        dataPlaneAppRetryRule);

    addStep(
        new PollCromwellRunSetStatusStep(
            flightBeanBag.getCbasService(),
            flightBeanBag.getSamService(),
            flightBeanBag.getImputationConfiguration()),
        dataPlaneAppRetryRule);
  }
}
