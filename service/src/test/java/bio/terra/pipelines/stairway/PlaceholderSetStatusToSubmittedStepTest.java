package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.JobStatusEnum;
import bio.terra.pipelines.service.JobsService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.time.Instant;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class PlaceholderSetStatusToSubmittedStepTest extends BaseContainerTest {

  @Autowired private JobsService jobsService;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void setStatus_doStep_success() throws InterruptedException {
    Instant beforeTimeSubmitted = Instant.now(); // for checking timeSubmitted
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    var setStatusStep = new PlaceholderSetStatusToSubmittedStep(jobsService);
    var result = setStatusStep.doStep(flightContext);

    Instant afterTimeSubmitted = Instant.now(); // for checking timeSubmitted

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // get info from the flight context to run checks
    FlightMap workingMap = flightContext.getWorkingMap();
    Instant timeSubmitted = workingMap.get(CreateJobFlightMapKeys.TIME_SUBMITTED, Instant.class);

    // make sure the status and time submitted were written to the working map
    assertEquals(
        JobStatusEnum.SUBMITTED.name(),
        workingMap.get(CreateJobFlightMapKeys.STATUS, String.class));
    // we can't check the exact time, but we can check that it's between the before and after times
    assertNotNull(timeSubmitted);
    assertTrue(beforeTimeSubmitted.isBefore(timeSubmitted));
    assertTrue(afterTimeSubmitted.isAfter(timeSubmitted));
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void setStatus_undoStep_success() throws InterruptedException {
    var setStatusStep = new PlaceholderSetStatusToSubmittedStep(jobsService);
    var result = setStatusStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
