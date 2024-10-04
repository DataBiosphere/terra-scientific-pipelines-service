package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.FlightUtils.inputParametersContainTrue;
import static bio.terra.pipelines.common.utils.FlightUtils.isContextInvalid;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
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

    if (inputParametersContainTrue(
            context.getInputParameters(),
            RunImputationJobFlightMapKeys.DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK)
        && context.getFlightStatus() != FlightStatus.SUCCESS) {
      logger.info(
          "Flight has status {}, setting PipelineRun status to FAILED", context.getFlightStatus());

      // set PipelineRun status to FAILED
      pipelineRunsService.markPipelineRunFailed(
          UUID.fromString(context.getFlightId()),
          context.getInputParameters().get(JobMapKeys.USER_ID.getKeyName(), String.class));
    }

    return HookAction.CONTINUE;
  }
}
