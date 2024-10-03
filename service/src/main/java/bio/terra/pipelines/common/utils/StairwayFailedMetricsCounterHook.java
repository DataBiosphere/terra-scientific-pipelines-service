package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.FlightUtils.getFlightClassFromString;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.imputation.PipelineRunTypeFlight;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.HookAction;
import bio.terra.stairway.StairwayHook;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

/**
 * A {@link StairwayHook} that logs a pipeline failure to Prometheus metrics upon flight failure.
 */
@Component
public class StairwayFailedMetricsCounterHook implements StairwayHook {
  private static final Logger logger =
      LoggerFactory.getLogger(StairwayFailedMetricsCounterHook.class);

  @Override
  public HookAction endFlight(FlightContext context) {
    if (isContextInvalid(context)) {
      return HookAction.CONTINUE;
    }

    Class<?> flightClass = getFlightClassFromString(context.getFlightClassName());
    if (flightClass == null) {
      logger.warn(
          "Failed to interpret flight class {}, skipping metrics hook",
          context.getFlightClassName());
      return HookAction.CONTINUE;
    } else if (flightClass.isInstance(PipelineRunTypeFlight.class)
        && context.getFlightStatus() != FlightStatus.SUCCESS) {
      logger.info("Flight failed, incrementing failed flight counter");

      // increment failed runs counter metric
      PipelinesEnum pipelinesEnum =
          PipelinesEnum.valueOf(
              context
                  .getInputParameters()
                  .get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
      MetricsUtils.incrementPipelineRunFailed(pipelinesEnum);
    }
    return HookAction.CONTINUE;
  }

  private boolean isContextInvalid(FlightContext context) {
    if (context == null || context.getWorkingMap() == null) {
      logger.warn("Flight context or working map null, skipping metrics hook");
      return true;
    }

    return false;
  }
}
