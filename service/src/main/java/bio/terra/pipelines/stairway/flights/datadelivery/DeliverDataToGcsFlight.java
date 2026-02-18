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
      new RetryRuleFixedInterval(/* intervalSeconds= */ 2, /* maxCount= */ 3);

  // addStep is protected in Flight, so make an override that is public
  @Override
  public void addStep(Step step, RetryRule retryRule) {
    super.addStep(step, retryRule);
  }

  public DeliverDataToGcsFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
    final FlightBeanBag flightBeanBag = FlightBeanBag.getFromObject(beanBag);

    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        JobMapKeys.PIPELINE_ID,
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH,
        DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    // Validate that the user has access to the destination GCS bucket
    // addStep(TODO TSPS-765)

    // Validate that the Teaspoons service account has access to the destination GCS bucket
    // addStep(TODO TSPS-765)

    // Copy the outputs to the user's specified destination
    addStep(
        new DeliverOutputFilesToGcsStep(
            flightBeanBag.getPipelineRunsService(),
            flightBeanBag.getPipelineInputsOutputsService()),
        gcsRetryRule);

    // Delete the outputs from the source workspace bucket
    // addStep(TODO);

    // Mark the pipelineRun as having completed data delivery, so outputs are no longer downloadable
    // and we have a record of where they went
    // addStep(TODO)
  }
}
