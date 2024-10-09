package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.FlightUtils.inputParametersContainTrue;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.HookAction;
import bio.terra.stairway.StairwayHook;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

/**
 * A {@link StairwayHook} that updates the PipelineRun status to FAILED on flight failure.
 *
 * <p>This the endFlight hook action will only run if the flight's input parameters contain the
 * JobMapKeys key for DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK and the flight's status is not SUCCESS.
 *
 * <p>The JobMapKeys key for USER_ID is required to set the PipelineRun status to FAILED.
 */
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

    if (inputParametersContainTrue(
            context.getInputParameters(), JobMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK)
        && context.getFlightStatus() != FlightStatus.SUCCESS) {
      logger.info(
          "Flight has status {}, setting PipelineRun status to FAILED", context.getFlightStatus());

      // set PipelineRun status to FAILED
      pipelineRunsService.markPipelineRunFailed(
          UUID.fromString(context.getFlightId()),
          context.getInputParameters().get(JobMapKeys.USER_ID, String.class));
    }

    return HookAction.CONTINUE;
  }
}
