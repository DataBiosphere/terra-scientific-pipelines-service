package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.service.JobsService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.time.Instant;
import java.util.Objects;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.RetryException;

public class WriteJobToDbStep implements Step {
  private final JobsService jobsService;
  private final Logger logger = LoggerFactory.getLogger(WriteJobToDbStep.class);

  public WriteJobToDbStep(JobsService jobsService) {
    this.jobsService = jobsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        CreateJobFlightMapKeys.PIPELINE_ID,
        CreateJobFlightMapKeys.PIPELINE_VERSION,
        CreateJobFlightMapKeys.SUBMITTING_USER_ID);

    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, CreateJobFlightMapKeys.STATUS);

    UUID writtenJobUUID =
        jobsService.writeJobToDb(
            UUID.fromString(flightContext.getFlightId()),
            inputParameters.get(CreateJobFlightMapKeys.SUBMITTING_USER_ID, String.class),
            inputParameters.get(CreateJobFlightMapKeys.PIPELINE_ID, String.class),
            inputParameters.get(CreateJobFlightMapKeys.PIPELINE_VERSION, String.class),
            workingMap.get(CreateJobFlightMapKeys.TIME_SUBMITTED, Instant.class),
            workingMap.get(CreateJobFlightMapKeys.STATUS, String.class),
            Objects.requireNonNull(
                inputParameters.get(CreateJobFlightMapKeys.PIPELINE_INPUTS, Object.class)));

    logger.info("Wrote job to db with id: {}", writtenJobUUID);
    workingMap.put(StairwayJobMapKeys.RESPONSE.getKeyName(), writtenJobUUID);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // TODO what do we need to do here? delete the row?
    return StepResult.getStepResultSuccess();
  }
}
