package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightDebugInfo;
import java.time.Duration;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class StairwayJobServiceTest extends BaseContainerTest {

  @Autowired StairwayJobService stairwayJobService;

  /**
   * Reset the {@link StairwayJobService} {@link FlightDebugInfo} after each test so that future
   * submissions aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    stairwayJobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void emptyStringJobIdTest() {
    String testJobId = "";
    assertThrows(
        InvalidStairwayJobIdException.class,
        () ->
            stairwayJobService
                .newJob()
                .description("description")
                .jobId(testJobId)
                .flightClass(StairwayJobServiceTestFlight.class));
  }

  @Test
  void whitespaceStringJobIdTest() {
    String testJobId = "\t ";

    assertThrows(
        InvalidStairwayJobIdException.class,
        () ->
            stairwayJobService
                .newJob()
                .description("description")
                .jobId(testJobId)
                .flightClass(StairwayJobServiceTestFlight.class));
  }

  @Test
  void submit_duplicateFlightId() {
    StairwayJobBuilder duplicateJobToSubmit =
        stairwayJobService
            .newJob()
            .description("description")
            .jobId("duplicateFlightId")
            .flightClass(StairwayJobServiceTestFlight.class);

    duplicateJobToSubmit.submit();

    assertThrows(DuplicateStairwayJobIdException.class, () -> duplicateJobToSubmit.submit());
  }

  @Test
  void testBadIdRetrieveJob() {
    assertThrows(
        StairwayJobNotFoundException.class, () -> stairwayJobService.retrieveJob("abcdef"));
  }

  @Test
  void testBadIdRetrieveResult() {
    assertThrows(
        StairwayJobNotFoundException.class,
        () -> stairwayJobService.retrieveJobResult("abcdef", Object.class));
  }

  @Test
  void setFlightDebugInfoForTest() throws InterruptedException {
    // Set a FlightDebugInfo so that any job submission should fail on the last step.
    stairwayJobService.setFlightDebugInfoForTest(
        FlightDebugInfo.newBuilder().lastStepFailure(true).build());

    String jobId = runFlight("fail for FlightDebugInfo");
    assertThrows(
        InvalidResultStateException.class,
        () -> stairwayJobService.retrieveJobResult(jobId, String.class));
  }

  // Submit a flight; wait for it to finish; return the flight id
  // Use the jobId defaulting in the JobBuilder
  private String runFlight(String description) throws InterruptedException {
    String jobId =
        stairwayJobService
            .newJob()
            .description(description)
            .flightClass(StairwayJobServiceTestFlight.class)
            .submit();
    StairwayTestUtils.pollUntilComplete(
        jobId, stairwayJobService.getStairway(), Duration.ofSeconds(1), Duration.ofSeconds(10));
    return jobId;
  }

  private String makeDescription(int ii) {
    return String.format("flight%d", ii);
  }
}
