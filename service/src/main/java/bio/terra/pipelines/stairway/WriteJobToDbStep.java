package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.service.ImputationService;
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
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        StairwayJobMapKeys.USER_ID.getKeyName(),
        StairwayJobMapKeys.PIPELINE_ID.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_VERSION);

    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationJobFlightMapKeys.STATUS);

    UUID writtenJobUUID =
        imputationService.writeJobToDb(
            UUID.fromString(flightContext.getFlightId()),
            inputParameters.get(StairwayJobMapKeys.USER_ID.getKeyName(), String.class),
            inputParameters.get(RunImputationJobFlightMapKeys.PIPELINE_VERSION, String.class),
            Objects.requireNonNull(
                inputParameters.get(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, Object.class)));

    logger.info("Wrote job to db with id: {}", writtenJobUUID);
    workingMap.put(StairwayJobMapKeys.RESPONSE.getKeyName(), writtenJobUUID);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // nothing to undo; keep the job in the database
    return StepResult.getStepResultSuccess();
  }
}
