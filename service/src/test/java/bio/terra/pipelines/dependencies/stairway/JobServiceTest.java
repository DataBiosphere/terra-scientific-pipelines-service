package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightDebugInfo;
import bio.terra.stairway.FlightState;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobServiceTest extends BaseEmbeddedDbTest {

  @Autowired JobService jobService;

  private static final PipelinesEnum PIPELINE_NAME = PipelinesEnum.IMPUTATION_BEAGLE;
  public static final Long PIPELINE_ID = 1L;
  private static final String TEST_USER_ID = TestUtils.TEST_USER_ID_1;
  private static final UUID TEST_JOB_UUID = TestUtils.TEST_NEW_UUID;

  /**
   * Reset the {@link JobService} {@link FlightDebugInfo} after each test so that future submissions
   * aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    jobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void submitDuplicateFlightId() throws InterruptedException {
    UUID jobId = TEST_JOB_UUID;

    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, "job for submit_duplicateFlightId() test")
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    jobToSubmit.submit();

    JobBuilder duplicateJob =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, "second job for submit_duplicateFlightId() test")
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    StairwayTestUtils.pollUntilComplete(jobId, jobService.getStairway(), 10L);

    assertThrows(DuplicateJobIdException.class, duplicateJob::submit);
  }

  @Test
  void submitSuccess() {
    // this tests the setters in JobBuilder
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(TEST_JOB_UUID)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, "job for submit_success() test")
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    // calling submit will run populateInputParameters() and validateRequiredInputs()
    assertDoesNotThrow(jobToSubmit::submit);
  }

  @Test
  void submitMissingFlightClass() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(TEST_JOB_UUID)
            .addParameter(
                JobMapKeys.DESCRIPTION, "description for submit_missingFlightClass() test")
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field for flight construction: jobId");
  }

  @Test
  void submitMissingUserId() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(TEST_JOB_UUID)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, "description for submit_missingUserId() test")
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field(s) for flight construction: userId");
  }

  @Test
  void submitNullRequiredField() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(TEST_JOB_UUID)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.USER_ID, null)
            .addParameter(JobMapKeys.DESCRIPTION, "description for submit_nullRequiredField() test")
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME);

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field(s) for flight construction: userId");
  }

  @Test
  void submitBlankRequiredField() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(TEST_JOB_UUID)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.USER_ID, "")
            .addParameter(JobMapKeys.DESCRIPTION, "description for submit_nullRequiredField() test")
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field(s) for flight construction: userId");
  }

  @Test
  void submitMissingPipelineId() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(TEST_JOB_UUID)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, "description for submit_missingPipelineId() test")
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME);

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field(s) for flight construction: pipelineId");
  }

  @Test
  void submitMissingJobId() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, "description for submit_missingJobId() test")
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field for flight construction: jobId");
  }

  @Test
  void submitMissingMultipleFields() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION, "description for submit_missingMultipleFields() test");

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field(s) for flight construction: userId, pipelineId");
  }

  @Test
  void submitMissingDescriptionOk() {
    UUID jobId = TEST_JOB_UUID;
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.USER_ID, TEST_USER_ID)
            .addParameter(JobMapKeys.PIPELINE_NAME, PIPELINE_NAME)
            .addParameter(JobMapKeys.PIPELINE_ID, PIPELINE_ID);

    // calling submit will run populateInputParameters() and validateRequiredInputs()
    assertDoesNotThrow(jobToSubmit::submit);
  }

  @Test
  void retrieveJobBadId() {
    UUID jobId = TEST_JOB_UUID;
    assertThrows(
        JobNotFoundException.class, () -> jobService.retrieveJob(jobId, TEST_USER_ID, null));
  }

  @Test
  void retrieveJobResultBadId() {
    UUID jobId = TEST_JOB_UUID;
    assertThrows(
        JobNotFoundException.class, () -> jobService.retrieveJobResult(jobId, Object.class));
  }

  /* Note: we currently only have one pipeline: Imputation. when we add the next pipeline,
  we should update this test with some instances of that pipeline as well. */
  @Test
  void enumerateJobsPipelineIdImputation() throws InterruptedException {
    // create two Imputation jobs
    UUID firstJobId = UUID.randomUUID();
    UUID secondJobId = UUID.randomUUID();
    // tests
    runFlight(firstJobId, TEST_USER_ID, PIPELINE_NAME, PIPELINE_ID, "imputation flight 1");
    runFlight(secondJobId, TEST_USER_ID, PIPELINE_NAME, PIPELINE_ID, "imputation flight 2");

    EnumeratedJobs jobs = jobService.enumerateJobs(TEST_USER_ID, 10, null, PIPELINE_NAME);
    assertEquals(2, jobs.getTotalResults());
  }

  @Test
  void enumerateJobsCorrectUserIsolation() throws InterruptedException {
    // create a job for the first user and verify that it shows up
    runFlight(TEST_JOB_UUID, TEST_USER_ID, PIPELINE_NAME, PIPELINE_ID, "first user's flight");
    EnumeratedJobs jobsUserOne = jobService.enumerateJobs(TEST_USER_ID, 10, null, null);
    assertEquals(1, jobsUserOne.getTotalResults());

    // create a job for the second user
    UUID jobIdUserTwo = UUID.randomUUID();
    String testUserId2 = TestUtils.TEST_USER_ID_2;
    runFlight(jobIdUserTwo, testUserId2, PIPELINE_NAME, PIPELINE_ID, "second user's flight");

    // Verify that the old userid still shows only 1 record
    EnumeratedJobs jobsUserOneAgain = jobService.enumerateJobs(TEST_USER_ID, 10, null, null);
    assertEquals(1, jobsUserOneAgain.getTotalResults());

    // Verify the new user's id shows a single job as well
    EnumeratedJobs jobsUserTwo = jobService.enumerateJobs(testUserId2, 10, null, null);
    assertEquals(1, jobsUserTwo.getTotalResults());
  }

  @Test
  void retrieveJobCorrectUserIsolation() throws InterruptedException {
    // create a job for the first user and verify that it shows up
    UUID jobIdUser1 = TEST_JOB_UUID;
    runFlight(jobIdUser1, TEST_USER_ID, PIPELINE_NAME, PIPELINE_ID, "first user's flight");
    FlightState user1job = jobService.retrieveJob(jobIdUser1, TEST_USER_ID, null);
    assertEquals(jobIdUser1.toString(), user1job.getFlightId());

    // make sure that user 2 doesn't have access to user 1's job
    assertThrows(
        JobUnauthorizedException.class,
        () -> jobService.retrieveJob(jobIdUser1, TestUtils.TEST_USER_ID_2, null));
  }

  @Test
  void retrieveJobWithPipelineName() throws InterruptedException {
    // create an imputation job for the first user and verify that it shows up
    UUID jobIdUser1 = TEST_JOB_UUID;
    runFlight(jobIdUser1, TEST_USER_ID, PIPELINE_NAME, PIPELINE_ID, "first user's flight");
    FlightState user1job = jobService.retrieveJob(jobIdUser1, TEST_USER_ID, PIPELINE_NAME);
    assertEquals(jobIdUser1.toString(), user1job.getFlightId());
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
    return runFlight(UUID.randomUUID(), TEST_USER_ID, PIPELINE_NAME, PIPELINE_ID, description);
  }

  // Submit a flight; wait for it to finish; return the flight id
  private UUID runFlight(
      UUID jobId, String userId, PipelinesEnum pipelineName, Long pipelineId, String description)
      throws InterruptedException {
    UUID submittedJobId =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION, description)
            .addParameter(JobMapKeys.USER_ID, userId)
            .addParameter(JobMapKeys.PIPELINE_NAME, pipelineName)
            .addParameter(JobMapKeys.PIPELINE_ID, pipelineId)
            .submit();
    StairwayTestUtils.pollUntilComplete(submittedJobId, jobService.getStairway(), 10L);
    return submittedJobId;
  }
}
