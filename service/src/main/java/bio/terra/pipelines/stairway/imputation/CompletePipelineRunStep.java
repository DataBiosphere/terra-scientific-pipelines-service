package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CompletePipelineRunStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final Logger logger = LoggerFactory.getLogger(CompletePipelineRunStep.class);

  public CompletePipelineRunStep(PipelineRunsService pipelineRunsService) {
    this.pipelineRunsService = pipelineRunsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID.getKeyName());

    UUID jobId = UUID.fromString(flightContext.getFlightId());
    String userId = inputParameters.get(JobMapKeys.USER_ID.getKeyName(), String.class);

    pipelineRunsService.markPipelineRunSuccess(jobId, userId);

    logger.info("Marked run {} as a success", jobId);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
