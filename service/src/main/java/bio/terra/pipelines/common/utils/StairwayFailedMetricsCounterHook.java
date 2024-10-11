package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.FlightUtils.flightMapKeyIsTrue;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.HookAction;
import bio.terra.stairway.StairwayHook;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

/**
 * A {@link StairwayHook} that logs a pipeline failure to Prometheus metrics upon flight failure.
 *
 * <p>This hook action will only run if the flight's input parameters contain the JobMapKeys key for
 * DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK and the flight's status is not SUCCESS.
 *
 * <p>The JobMapKeys key for PIPELINE_NAME is required to increment the failed flight.
 */
@Component
public class StairwayFailedMetricsCounterHook implements StairwayHook {
  private static final Logger logger =
      LoggerFactory.getLogger(StairwayFailedMetricsCounterHook.class);

  @Override
  public HookAction endFlight(FlightContext context) {

    if (flightMapKeyIsTrue(
            context.getInputParameters(), JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK)
        && context.getFlightStatus() != FlightStatus.SUCCESS) {
      logger.info(
          "Flight has status {}, incrementing failed flight counter", context.getFlightStatus());

      // increment failed runs counter metric
      PipelinesEnum pipelinesEnum =
          PipelinesEnum.valueOf(
              context.getInputParameters().get(JobMapKeys.PIPELINE_NAME, String.class));
      MetricsUtils.incrementPipelineRunFailed(pipelinesEnum);
    }
    return HookAction.CONTINUE;
  }
}
