package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
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
        RunImputationGcpJobFlightMapKeys.PIPELINE_ID,
        RunImputationGcpJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
        RunImputationGcpJobFlightMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
        RunImputationGcpJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        RunImputationGcpJobFlightMapKeys.CONTROL_WORKSPACE_PROJECT,
        RunImputationGcpJobFlightMapKeys.CONTROL_WORKSPACE_NAME,
        RunImputationGcpJobFlightMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_URL,
        RunImputationGcpJobFlightMapKeys.WDL_METHOD_NAME,
        JobMapKeys.RESULT_PATH.getKeyName());

    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    MetricsUtils.incrementPipelineRun(pipelinesEnum);

    addStep(
        new PrepareImputationInputsStep(
            flightBeanBag.getPipelinesService(), flightBeanBag.getImputationConfiguration()),
        dbRetryRule);

    addStep(new CompletePipelineRunStep(flightBeanBag.getPipelineRunsService()), dbRetryRule);
  }
}