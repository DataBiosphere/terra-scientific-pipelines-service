package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightDebugInfo;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class StairwayJobServiceTest extends BaseContainerTest {

  @Autowired StairwayJobService stairwayJobService;

  private static final String testPipelineId = MockMvcUtils.TEST_PIPELINE_ID_1;
  private static final String testDescription = "description";
  private static final String testRequest = "request";

  private static final String testPipelineVersion = MockMvcUtils.TEST_PIPELINE_VERSION_1;
  private static final String testUserId = MockMvcUtils.TEST_USER_ID_1;

  private final Object testPipelineInputs = MockMvcUtils.TEST_PIPELINE_INPUTS;

  /**
   * Reset the {@link StairwayJobService} {@link FlightDebugInfo} after each test so that future
   * submissions aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    stairwayJobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void submit_duplicateFlightId() throws InterruptedException {
    UUID jobId = MockMvcUtils.TEST_EXISTING_UUID;

    StairwayJobBuilder jobToSubmit =
        stairwayJobService
            .newJob()
            .description(testDescription)
            .request(testRequest)
            .pipelineId(testPipelineId)
            .jobId(jobId)
            .flightClass(StairwayJobServiceTestFlight.class);

    jobToSubmit.submit();

    StairwayJobBuilder duplicateJob =
        stairwayJobService.newJob().jobId(jobId).flightClass(StairwayJobServiceTestFlight.class);

    StairwayTestUtils.pollUntilComplete(jobId, stairwayJobService.getStairway(), 10L);

    assertThrows(DuplicateStairwayJobIdException.class, duplicateJob::submit);
  }

  @Test
  void submit_success() {
    // this tests the setters in StairwayJobBuilder
    StairwayJobBuilder jobToSubmit =
        stairwayJobService
            .newJob()
            .jobId(MockMvcUtils.TEST_NEW_UUID)
            .flightClass(StairwayJobServiceTestFlight.class)
            .description(testDescription)
            .request(testRequest)
            .pipelineId(testPipelineId)
            .pipelineVersion(testPipelineVersion)
            .submittingUserId(testUserId)
            .pipelineInputs(testPipelineInputs);

    // calling submit will run populateInputParameters()
    assertDoesNotThrow(jobToSubmit::submit);
  }

  @Test
  void submit_missingFlightClass() {
    StairwayJobBuilder jobToSubmit =
        stairwayJobService
            .newJob()
            .jobId(MockMvcUtils.TEST_NEW_UUID)
            .description(testDescription)
            .request(testRequest)
            .pipelineId(testPipelineId);

    assertThrows(MissingRequiredFieldException.class, jobToSubmit::submit);
  }

  @Test
  void retrieveJob_badId() {
    UUID unsubmittedJobId = UUID.randomUUID();
    assertThrows(
        StairwayJobNotFoundException.class, () -> stairwayJobService.retrieveJob(unsubmittedJobId));
  }

  @Test
  void retrieveJobResult_badId() {
    UUID unsubmittedJobId = UUID.randomUUID();
    assertThrows(
        StairwayJobNotFoundException.class,
        () -> stairwayJobService.retrieveJobResult(unsubmittedJobId, Object.class));
  }

  @Test
  void setFlightDebugInfoForTest() throws InterruptedException {
    // Set a FlightDebugInfo so that any job submission should fail on the last step.
    stairwayJobService.setFlightDebugInfoForTest(
        FlightDebugInfo.newBuilder().lastStepFailure(true).build());

    UUID jobId = runFlight("fail for FlightDebugInfo");
    assertThrows(
        InvalidResultStateException.class,
        () -> stairwayJobService.retrieveJobResult(jobId, UUID.class));
  }

  // Submit a flight; wait for it to finish; return the flight id
  // Use the jobId defaulting in the JobBuilder
  private UUID runFlight(String description) throws InterruptedException {
    UUID jobId =
        stairwayJobService
            .newJob()
            .jobId(UUID.randomUUID())
            .description(description)
            .flightClass(StairwayJobServiceTestFlight.class)
            .submit();
    StairwayTestUtils.pollUntilComplete(jobId, stairwayJobService.getStairway(), 10L);
    return jobId;
  }
}
