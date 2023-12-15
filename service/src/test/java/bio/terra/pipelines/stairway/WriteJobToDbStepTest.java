package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.CommonJobStatusEnum;
import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.service.JobsService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.time.Instant;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class WriteJobToDbStepTest extends BaseContainerTest {

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
  void writeJob_doStep_success() throws InterruptedException {
    // setup
    String testJobId = UUID.randomUUID().toString();
    when(flightContext.getFlightId()).thenReturn(testJobId);

    Instant timeSubmitted = Instant.now();
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    flightContext
        .getWorkingMap()
        .put(CreateJobFlightMapKeys.STATUS, CommonJobStatusEnum.SUBMITTED.name());
    flightContext.getWorkingMap().put(CreateJobFlightMapKeys.TIME_SUBMITTED, timeSubmitted);

    // do the step
    var writeJobStep = new WriteJobToDbStep(jobsService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();
    FlightMap workingMap = flightContext.getWorkingMap();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(testJobId, workingMap.get(StairwayJobMapKeys.RESPONSE.getKeyName(), String.class));

    // make sure the job was written to the db
    Job writtenJob =
        jobsService.getJob(
            inputParams.get(CreateJobFlightMapKeys.SUBMITTING_USER_ID, String.class),
            inputParams.get(CreateJobFlightMapKeys.PIPELINE_ID, String.class),
            UUID.fromString(testJobId));
    assertEquals(MockMvcUtils.TEST_PIPELINE_VERSION_1, writtenJob.getPipelineVersion());
    assertEquals(CommonJobStatusEnum.SUBMITTED.name(), writtenJob.getStatus());
    assertEquals(timeSubmitted, writtenJob.getTimeSubmitted());
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void writeJob_undoStep_success() throws InterruptedException {
    var writeJobStep = new WriteJobToDbStep(jobsService);
    var result = writeJobStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
