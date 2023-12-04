package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidStairwayJobIdException;
import bio.terra.pipelines.dependencies.stairway.exception.StairwayJobNotFoundException;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightDebugInfo;
import java.time.Duration;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class StairwayJobServiceTest extends BaseContainerTest {

  @Autowired private StairwayJobService stairwayJobService;
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

  //  // Resets the application context before retrieveTest to make sure that the job service does
  // not
  //  // have some failed jobs left over from other tests.
  //  @DirtiesContext(methodMode = DirtiesContext.MethodMode.BEFORE_METHOD)
  //  @Test
  //  void retrieveTest() {
  //    // We perform 7 flights and then retrieve and enumerate them.
  //    // The fids list should be in exactly the same order as the database ordered by submit time.
  //
  //    List<String> jobIds1 = new ArrayList<>();
  //    UUID workspace1 = WorkspaceUnitTestUtils.createWorkspaceWithoutCloudContext(workspaceDao);
  //    for (int i = 0; i < 3; i++) {
  //      String jobId = runFlight(workspace1, makeDescription(i));
  //      jobIds1.add(jobId);
  //    }
  //
  //    List<String> jobIds2 = new ArrayList<>();
  //    UUID workspace2 = WorkspaceUnitTestUtils.createWorkspaceWithoutCloudContext(workspaceDao);
  //
  //    for (int i = 0; i < 4; i++) {
  //      String jobId = runFlight(workspace2, makeDescription(i));
  //      jobIds2.add(jobId);
  //    }
  //
  //    // Test single retrieval
  //    testSingleRetrieval(jobIds1);
  //
  //    // Test result retrieval - the body should be the description string
  //    testResultRetrieval(jobIds1);
  //
  //    // Retrieve each workspace
  //    testEnumCount(jobIds1, workspace1, 0, 3, 100, null);
  //    testEnumCount(jobIds2, workspace2, 0, 4, 100, null);
  //
  //    // Test page token
  //    String pageToken = testEnumCount(jobIds2, workspace2, 0, 2, 2, null);
  //    pageToken = testEnumCount(jobIds2, workspace2, 2, 2, 2, pageToken);
  //    testEnumCount(jobIds2, workspace2, 4, 0, 2, pageToken);
  //  }

  //  private void testSingleRetrieval(List<String> fids) {
  //    FlightState flightState = stairwayJobService.retrieveJob(fids.get(2));
  //    ApiJobReport response = jobApiUtils.mapFlightStateToApiJobReport(flightState);
  //    assertThat(response, notNullValue());
  //    validateJobReport(response, 2, fids);
  //  }
  //
  //  private void testResultRetrieval(List<String> fids) {
  //    stairwayJobService.JobResultOrException<String> resultHolder =
  //        stairwayJobService.retrieveJobResult(fids.get(2), String.class);
  //
  //    assertNull(resultHolder.getException());
  //    assertEquals(resultHolder.getResult(), makeDescription(2));
  //  }
  //
  //  // Enumerate and make sure we got the number we expected
  //  // Validate the result is what we expect
  //  private String testEnumCount(
  //      List<String> fids,
  //      UUID workspaceId,
  //      int expectedOffset,
  //      int expectedCount,
  //      int limit,
  //      String pageToken) {
  //    EnumeratedJobs jobList =
  //        stairwayJobService.enumerateJobs(workspaceId, limit, pageToken, null, null, null, null);
  //    assertNotNull(jobList);
  //    assertEquals(expectedCount, jobList.getResults().size());
  //    int index = expectedOffset;
  //    for (EnumeratedJob job : jobList.getResults()) {
  //      ApiJobReport jobReport = jobApiUtils.mapFlightStateToApiJobReport(job.getFlightState());
  //      validateJobReport(jobReport, index, fids);
  //      index++;
  //    }
  //
  //    return jobList.getPageToken();
  //  }

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
