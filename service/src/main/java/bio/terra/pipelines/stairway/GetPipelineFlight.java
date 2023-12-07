package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.stairway.*;

/**
 * Flight for creation of a controlled resource. Some steps are resource-type-agnostic, and others
 * depend on the resource type. The latter must be passed in via the input parameters map with keys
 */
public class GetPipelineFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 1, /* maxCount= */ 5);

  // addStep is protected in Flight, so make an override that is public
  @Override
  public void addStep(Step step, RetryRule retryRule) {
    super.addStep(step, retryRule);
  }

  public GetPipelineFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
    final FlightBeanBag flightBeanBag = FlightBeanBag.getFromObject(beanBag);

    FlightUtils.validateRequiredEntries(inputParameters, GetPipelineFlightMapKeys.PIPELINE_ID);

    // query the Pipelines table
    addStep(new GetPipelineStep(flightBeanBag.getPipelinesService()), dbRetryRule);
  }
}
