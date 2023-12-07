package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightDebugInfo;
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
  void newJob_emptyStringJobId() {
    String testJobId = "";

    StairwayJobBuilder baseJob = stairwayJobService.newJob();

    assertThrows(InvalidStairwayJobIdException.class, () -> baseJob.jobId(testJobId));
  }

  @Test
  void newJob_whitespaceJobId() {
    String testJobId = "\t ";

    StairwayJobBuilder baseJob = stairwayJobService.newJob();

    assertThrows(InvalidStairwayJobIdException.class, () -> baseJob.jobId(testJobId));
  }

  @Test
  void submit_duplicateFlightId() throws InterruptedException {
    String jobId = "duplicateFlightId";

    StairwayJobBuilder jobToSubmit =
        stairwayJobService
            .newJob()
            .description("description")
            .request("request")
            .pipelineId("pipelineId")
            .jobId(jobId)
            .flightClass(StairwayJobServiceTestFlight.class);

    jobToSubmit.submit();

    StairwayJobBuilder duplicateJob =
        stairwayJobService.newJob().jobId(jobId).flightClass(StairwayJobServiceTestFlight.class);

    StairwayTestUtils.pollUntilComplete(jobId, stairwayJobService.getStairway(), 10L);

    assertThrows(DuplicateStairwayJobIdException.class, () -> duplicateJob.submit());
  }

  @Test
  void submit_missingFlightClass() {
    StairwayJobBuilder jobToSubmit =
        stairwayJobService
            .newJob()
            .description("description")
            .request("request")
            .pipelineId("pipelineId");

    assertThrows(MissingRequiredFieldException.class, () -> jobToSubmit.submit());
  }

  @Test
  void retrieveJob_badId() {
    assertThrows(
        StairwayJobNotFoundException.class, () -> stairwayJobService.retrieveJob("abcdef"));
  }

  @Test
  void retrieveJobResult_badId() {
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
    StairwayTestUtils.pollUntilComplete(jobId, stairwayJobService.getStairway(), 10L);
    return jobId;
  }
}
