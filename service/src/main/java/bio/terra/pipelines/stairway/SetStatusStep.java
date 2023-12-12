package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.JobStatusEnum;
import bio.terra.pipelines.service.JobsService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.time.Instant;

/* This is a placeholder step to set the status in the working map; it will be replaced with real steps in future PRs */
public class SetStatusStep implements Step {
  private final JobsService jobsService;

  public SetStatusStep(JobsService jobsService) {
    this.jobsService = jobsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) throws InterruptedException {
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(CreateJobFlightMapKeys.STATUS, JobStatusEnum.SUBMITTED.name());

    Instant timeSubmitted = jobsService.getCurrentTimestamp();
    workingMap.put(CreateJobFlightMapKeys.TIME_SUBMITTED, timeSubmitted);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    return StepResult.getStepResultSuccess();
  }
}
