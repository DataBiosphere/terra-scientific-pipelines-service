package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.imputation.gcp.AddDataTableRowStep;
import bio.terra.pipelines.stairway.imputation.gcp.PollCromwellSubmissionStatusStep;
import bio.terra.pipelines.stairway.imputation.gcp.SubmitCromwellRunSetStep;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.RetryRule;
import bio.terra.stairway.RetryRuleFixedInterval;
import bio.terra.stairway.Step;

public class RunImputationGcpJobFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 1, /* maxCount= */ 5);

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
        JobMapKeys.USER_ID.getKeyName(),
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_ID,
        RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
        RunImputationJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
        RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_PROJECT,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_URL,
        RunImputationJobFlightMapKeys.WDL_METHOD_NAME,
        JobMapKeys.RESULT_PATH.getKeyName());

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    addStep(
        new PrepareImputationInputsStep(
            flightBeanBag.getPipelinesService(), flightBeanBag.getImputationConfiguration()),
        dbRetryRule);

    addStep(
        new AddDataTableRowStep(flightBeanBag.getRawlsService(), flightBeanBag.getSamService()));

    addStep(
        new SubmitCromwellRunSetStep(
            flightBeanBag.getRawlsService(), flightBeanBag.getSamService()),
        dbRetryRule);

    addStep(
        new PollCromwellSubmissionStatusStep(
            flightBeanBag.getSamService(),
            flightBeanBag.getRawlsService(),
            flightBeanBag.getImputationConfiguration()),
        dbRetryRule);
  }
}
