package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.JobStatusEnum;
import bio.terra.pipelines.service.JobsService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.time.Instant;

public class PlaceholderSetStatusToSubmittedStep implements Step {
  private final JobsService jobsService;

  @SuppressWarnings(
      "java:S125") // this comment block will be removed once this is converted to a real step
  /* This is a placeholder step that only sets the status in the working map;
  it will be replaced with real steps in future PRs */
  public PlaceholderSetStatusToSubmittedStep(JobsService jobsService) {
    this.jobsService = jobsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) throws InterruptedException {

    // to add later: submit the workflow to CBAS

    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(CreateJobFlightMapKeys.STATUS, JobStatusEnum.SUBMITTED.name());

    Instant timeSubmitted = jobsService.getCurrentTimestamp();
    workingMap.put(CreateJobFlightMapKeys.TIME_SUBMITTED, timeSubmitted);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // nothing to undo yet
    return StepResult.getStepResultSuccess();
  }
}
