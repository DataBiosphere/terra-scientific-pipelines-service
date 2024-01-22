package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightDebugInfo;
import bio.terra.stairway.FlightState;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobServiceTest extends BaseContainerTest {

  @Autowired JobService jobService;

  private static final PipelinesEnum imputationPipelineId = PipelinesEnum.IMPUTATION;
  private static final String testUserId = TestUtils.TEST_USER_ID_1;
  private static final UUID newJobId = TestUtils.TEST_NEW_UUID;

  /**
   * Reset the {@link JobService} {@link FlightDebugInfo} after each test so that future submissions
   * aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    jobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void submit_duplicateFlightId() throws InterruptedException {
    // TSPS-128 will make these tests truly independent, and should update tests here to use the
    // newJobId
    UUID jobId = UUID.randomUUID(); // newJobId;

    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(), "job for submit_duplicateFlightId() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), imputationPipelineId);

    jobToSubmit.submit();

    JobBuilder duplicateJob =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "second job for submit_duplicateFlightId() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), imputationPipelineId);

    StairwayTestUtils.pollUntilComplete(jobId, jobService.getStairway(), 10L);

    assertThrows(DuplicateJobIdException.class, duplicateJob::submit);
  }

  @Test
  void submit_success() {
    // this tests the setters in JobBuilder
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(UUID.randomUUID()) // newJobId
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), "job for submit_success() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), imputationPipelineId);

    // calling submit will run populateInputParameters()
    assertDoesNotThrow(jobToSubmit::submit);
  }

  @Test
  void submit_missingFlightClass() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(UUID.randomUUID()) // newJobId
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_missingFlightClass() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), imputationPipelineId);

    assertThrows(MissingRequiredFieldException.class, jobToSubmit::submit);
  }

  @Test
  void submit_missingUserId() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(), "description for submit_missingUserId() test")
            .addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), imputationPipelineId);

    assertThrows(MissingRequiredFieldException.class, jobToSubmit::submit);
  }

  @Test
  void submit_missingPipelineId() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_missingPipelineId() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId);

    assertThrows(MissingRequiredFieldException.class, jobToSubmit::submit);
  }

  @Test
  void retrieveJob_badId() {
    UUID jobId = UUID.randomUUID(); // newJobId
    assertThrows(JobNotFoundException.class, () -> jobService.retrieveJob(jobId, testUserId));
  }

  @Test
  void retrieveJobResult_badId() {
    UUID jobId = UUID.randomUUID(); // newJobId
    assertThrows(
        JobNotFoundException.class, () -> jobService.retrieveJobResult(jobId, Object.class));
  }

  /* Note: we currently only have one pipeline: Imputation. when we add the next pipeline,
  we should update this test with some instances of that pipeline as well. */
  @Test
  void testEnumerateJobsPipelineIdImputation() throws InterruptedException {
    // create two Imputation jobs
    UUID firstJobId = UUID.randomUUID();
    UUID secondJobId = UUID.randomUUID();
    String newTestUserId =
        "anotherUserId"; // use testUserId once we implement TSPS-128 for effectively independent
    // tests
    runFlight(firstJobId, newTestUserId, imputationPipelineId, "imputation flight 1");
    runFlight(secondJobId, newTestUserId, imputationPipelineId, "imputation flight 2");

    EnumeratedJobs jobs = jobService.enumerateJobs(newTestUserId, 10, null, imputationPipelineId);
    assertEquals(2, jobs.getTotalResults());
  }

  @Test
  void testEnumerateJobsCorrectUserIsolation() throws InterruptedException {
    // create a job for the first user and verify that it shows up
    runFlight(newJobId, testUserId, imputationPipelineId, "first user's flight");
    EnumeratedJobs jobsUserOne = jobService.enumerateJobs(testUserId, 10, null, null);
    assertEquals(1, jobsUserOne.getTotalResults());

    // create a job for the second user
    UUID jobIdUserTwo = UUID.randomUUID();
    String testUserId2 = TestUtils.TEST_USER_ID_2;
    runFlight(jobIdUserTwo, testUserId2, imputationPipelineId, "second user's flight");

    // Verify that the old userid still shows only 1 record
    EnumeratedJobs jobsUserOneAgain = jobService.enumerateJobs(testUserId, 10, null, null);
    assertEquals(1, jobsUserOneAgain.getTotalResults());

    // Verify the new user's id shows a single job as well
    EnumeratedJobs jobsUserTwo = jobService.enumerateJobs(testUserId2, 10, null, null);
    assertEquals(1, jobsUserTwo.getTotalResults());
  }

  @Test
  void testRetrieveJobCorrectUserIsolation() throws InterruptedException {
    // create a job for the first user and verify that it shows up
    UUID jobIdUser1 = UUID.randomUUID(); // newJobId
    runFlight(jobIdUser1, testUserId, imputationPipelineId, "first user's flight");
    FlightState user1job = jobService.retrieveJob(jobIdUser1, testUserId);
    assertEquals(jobIdUser1.toString(), user1job.getFlightId());

    // make sure that user 2 doesn't have access to user 1's job
    assertThrows(
        JobUnauthorizedException.class,
        () -> jobService.retrieveJob(jobIdUser1, TestUtils.TEST_USER_ID_2));
  }

  @Test
  void setFlightDebugInfoForTest() throws InterruptedException {
    // Set a FlightDebugInfo so that any job submission should fail on the last step.
    jobService.setFlightDebugInfoForTest(
        FlightDebugInfo.newBuilder().lastStepFailure(true).build());

    UUID jobId = runFlight("fail for FlightDebugInfo");
    assertThrows(
        InvalidResultStateException.class, () -> jobService.retrieveJobResult(jobId, UUID.class));
  }

  // Submit a flight; wait for it to finish; return the flight id
  // using randomly generated flightId and the test userId
  private UUID runFlight(String description) throws InterruptedException {
    return runFlight(UUID.randomUUID(), testUserId, imputationPipelineId, description);
  }

  // Submit a flight; wait for it to finish; return the flight id
  private UUID runFlight(UUID jobId, String userId, PipelinesEnum pipelineId, String description)
      throws InterruptedException {
    UUID submittedJobId =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), userId)
            .addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), pipelineId)
            .submit();
    StairwayTestUtils.pollUntilComplete(submittedJobId, jobService.getStairway(), 10L);
    return submittedJobId;
  }
}
