package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.stairway.*;

public class RunImputationJobFlight extends Flight {

  /** Retry for short database operations which may fail due to transaction conflicts. */
  private final RetryRule dbRetryRule =
      new RetryRuleFixedInterval(/*intervalSeconds= */ 1, /* maxCount= */ 5);

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
        JobMapKeys.PIPELINE_ID.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_VERSION,
        RunImputationJobFlightMapKeys.PIPELINE_INPUTS);

    // this currently just sets the status to SUBMITTED and puts the current time into the working
    // map
    addStep(new PlaceholderSetStatusToSubmittedStep());

    // write the job metadata to the Jobs table
    addStep(new WriteJobToDbStep(flightBeanBag.getImputationService()), dbRetryRule);
  }
}
