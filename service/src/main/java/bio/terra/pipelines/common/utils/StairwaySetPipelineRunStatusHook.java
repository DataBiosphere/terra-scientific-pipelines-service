package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.FlightUtils.getFlightClassFromString;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.imputation.PipelineRunTypeFlight;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.HookAction;
import bio.terra.stairway.StairwayHook;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

/** A {@link StairwayHook} that updates the PipelineRun status to FAILED on flight failure. */
@Component
public class StairwaySetPipelineRunStatusHook implements StairwayHook {
  private static final Logger logger =
      LoggerFactory.getLogger(StairwaySetPipelineRunStatusHook.class);
  private final PipelineRunsService pipelineRunsService;

  public StairwaySetPipelineRunStatusHook(PipelineRunsService pipelineRunsService) {
    this.pipelineRunsService = pipelineRunsService;
  }

  @Override
  public HookAction endFlight(FlightContext context) {
    if (isContextInvalid(context)) {
      return HookAction.CONTINUE;
    }

    Class<?> flightClass = getFlightClassFromString(context.getFlightClassName());
    if (flightClass == null) {
      logger.warn(
          "Failed to interpret flight class {}, skipping update status hook",
          context.getFlightClassName());
      return HookAction.CONTINUE;
    } else if (flightClass.isInstance(PipelineRunTypeFlight.class)
        && context.getFlightStatus() != FlightStatus.SUCCESS) {
      logger.info(
          "Flight has status {}, setting PipelineRun status to FAILED", context.getFlightStatus());

      // set PipelineRun status to FAILED
      var inputParameters = context.getInputParameters();
      FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID.getKeyName());
      pipelineRunsService.markPipelineRunFailed(
          UUID.fromString(context.getFlightId()),
          inputParameters.get(JobMapKeys.USER_ID.getKeyName(), String.class));
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
