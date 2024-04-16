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

  private static final PipelinesEnum imputationPipelineName = PipelinesEnum.IMPUTATION_BEAGLE;
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
  void submitDuplicateFlightId() throws InterruptedException {
    UUID jobId = newJobId;

    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(), "job for submit_duplicateFlightId() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

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
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

    StairwayTestUtils.pollUntilComplete(jobId, jobService.getStairway(), 10L);

    assertThrows(DuplicateJobIdException.class, duplicateJob::submit);
  }

  @Test
  void submitSuccess() {
    // this tests the setters in JobBuilder
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), "job for submit_success() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

    // calling submit will run populateInputParameters() and validateRequiredInputs()
    assertDoesNotThrow(jobToSubmit::submit);
  }

  @Test
  void submitMissingFlightClass() {
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(newJobId)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_missingFlightClass() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

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
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(), "description for submit_missingUserId() test")
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

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
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), null)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_nullRequiredField() test")
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

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
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), "")
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_nullRequiredField() test")
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

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
            .jobId(newJobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_missingPipelineId() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId);

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
            .addParameter(
                JobMapKeys.DESCRIPTION.getKeyName(), "description for submit_missingJobId() test")
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);
    ;

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
                JobMapKeys.DESCRIPTION.getKeyName(),
                "description for submit_missingMultipleFields() test");
    ;

    assertThrows(
        MissingRequiredFieldException.class,
        jobToSubmit::submit,
        "Missing required field(s) for flight construction: userId, pipelineId");
  }

  @Test
  void submitMissingDescriptionOk() {
    UUID jobId = newJobId;
    JobBuilder jobToSubmit =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(JobServiceTestFlight.class)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName);

    // calling submit will run populateInputParameters() and validateRequiredInputs()
    assertDoesNotThrow(jobToSubmit::submit);
  }

  @Test
  void retrieveJobBadId() {
    UUID jobId = newJobId;
    assertThrows(JobNotFoundException.class, () -> jobService.retrieveJob(jobId, testUserId, null));
  }

  @Test
  void retrieveJobResultBadId() {
    UUID jobId = newJobId;
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
    runFlight(firstJobId, testUserId, imputationPipelineName, "imputation flight 1");
    runFlight(secondJobId, testUserId, imputationPipelineName, "imputation flight 2");

    EnumeratedJobs jobs = jobService.enumerateJobs(testUserId, 10, null, imputationPipelineName);
    assertEquals(2, jobs.getTotalResults());
  }

  @Test
  void enumerateJobsCorrectUserIsolation() throws InterruptedException {
    // create a job for the first user and verify that it shows up
    runFlight(newJobId, testUserId, imputationPipelineName, "first user's flight");
    EnumeratedJobs jobsUserOne = jobService.enumerateJobs(testUserId, 10, null, null);
    assertEquals(1, jobsUserOne.getTotalResults());

    // create a job for the second user
    UUID jobIdUserTwo = UUID.randomUUID();
    String testUserId2 = TestUtils.TEST_USER_ID_2;
    runFlight(jobIdUserTwo, testUserId2, imputationPipelineName, "second user's flight");

    // Verify that the old userid still shows only 1 record
    EnumeratedJobs jobsUserOneAgain = jobService.enumerateJobs(testUserId, 10, null, null);
    assertEquals(1, jobsUserOneAgain.getTotalResults());

    // Verify the new user's id shows a single job as well
    EnumeratedJobs jobsUserTwo = jobService.enumerateJobs(testUserId2, 10, null, null);
    assertEquals(1, jobsUserTwo.getTotalResults());
  }

  @Test
  void retrieveJobCorrectUserIsolation() throws InterruptedException {
    // create a job for the first user and verify that it shows up
    UUID jobIdUser1 = newJobId;
    runFlight(jobIdUser1, testUserId, imputationPipelineName, "first user's flight");
    FlightState user1job = jobService.retrieveJob(jobIdUser1, testUserId, null);
    assertEquals(jobIdUser1.toString(), user1job.getFlightId());

    // make sure that user 2 doesn't have access to user 1's job
    assertThrows(
        JobUnauthorizedException.class,
        () -> jobService.retrieveJob(jobIdUser1, TestUtils.TEST_USER_ID_2, null));
  }

  @Test
  void retrieveJobWithPipelineName() throws InterruptedException {
    // create an imputation job for the first user and verify that it shows up
    UUID jobIdUser1 = newJobId;
    runFlight(jobIdUser1, testUserId, imputationPipelineName, "first user's flight");
    FlightState user1job = jobService.retrieveJob(jobIdUser1, testUserId, imputationPipelineName);
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
    return runFlight(UUID.randomUUID(), testUserId, imputationPipelineName, description);
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
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), pipelineId)
            .submit();
    StairwayTestUtils.pollUntilComplete(submittedJobId, jobService.getStairway(), 10L);
    return submittedJobId;
  }
}
