package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.Objects;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.RetryException;

public class WriteJobToDbStep implements Step {
  private final ImputationService imputationService;
  private final Logger logger = LoggerFactory.getLogger(WriteJobToDbStep.class);

  public WriteJobToDbStep(ImputationService imputationService) {
    this.imputationService = imputationService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID.getKeyName(),
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_ID);

    UUID writtenJobUUID =
        imputationService.writeJobToDb(
            UUID.fromString(flightContext.getFlightId()),
            inputParameters.get(JobMapKeys.USER_ID.getKeyName(), String.class),
            inputParameters.get(RunImputationJobFlightMapKeys.PIPELINE_ID, Long.class),
            Objects.requireNonNull(
                inputParameters.get(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, Object.class)));

    logger.info("Wrote job to db with id: {}", writtenJobUUID);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // nothing to undo; keep the job in the database
    return StepResult.getStepResultSuccess();
  }
}
