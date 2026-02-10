package bio.terra.pipelines.stairway.flights.datadelivery;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.steps.datadelivery.DeliverOutputFilesToGcsStep;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.RetryRule;
import bio.terra.stairway.RetryRuleFixedInterval;
import bio.terra.stairway.Step;

public class DeliverDataToGcsFlight extends Flight {

  /** Retry rule for GCS operations which may fail due to transient issues. */
  private final RetryRule gcsRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 2, /* maxCount= */ 3);

  // addStep is protected in Flight, so make an override that is public
  @Override
  public void addStep(Step step, RetryRule retryRule) {
    super.addStep(step, retryRule);
  }

  public DeliverDataToGcsFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
    final FlightBeanBag flightBeanBag = FlightBeanBag.getFromObject(beanBag);

    // Validate required parameters for data delivery
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH,
        DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    // Add the single step to deliver output files to GCS
    addStep(new DeliverOutputFilesToGcsStep(flightBeanBag.getPipelineRunsService()), gcsRetryRule);
  }
}
