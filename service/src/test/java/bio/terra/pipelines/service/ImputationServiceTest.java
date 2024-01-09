package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.ImputationJob;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.ImputationJobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.time.Instant;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class ImputationServiceTest extends BaseContainerTest {

  @Autowired ImputationService imputationService;
  @Autowired ImputationJobsRepository imputationJobsRepository;
  @Autowired PipelineInputsRepository pipelineInputsRepository;

  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private ImputationJob createTestJobWithJobId(UUID jobId) {
    return createTestJobWithJobIdAndUser(jobId, testUserId);
  }

  private ImputationJob createTestJobWithJobIdAndUser(UUID jobId, String userId) {
    return new ImputationJob(jobId, userId, testPipelineVersion);
  }

  @Test
  void testGetTimestamp() {
    Instant timestamp1 = imputationService.getCurrentTimestamp();
    Instant timestamp2 = imputationService.getCurrentTimestamp();

    assertNotNull(timestamp1);
    assertNotNull(timestamp2);
    assertNotEquals(timestamp1, timestamp2);
    assertTrue(timestamp1.isBefore(timestamp2));
  }

  @Test
  void testCreateJobId() {
    UUID jobId1 = imputationService.createJobId();
    UUID jobId2 = imputationService.createJobId();

    assertNotNull(jobId1);
    assertNotNull(jobId2);
    assertNotEquals(jobId1, jobId2);
  }

  @Test
  void testWriteJobToDb() {
    List<ImputationJob> jobsDefault = imputationService.getImputationJobs(testUserId);
    // test data migration inserts one row by default
    assertEquals(1, jobsDefault.size());

    UUID savedUUID =
        imputationService.writeJobToDb(
            testJobId, testUserId, testPipelineVersion, testPipelineInputs);

    List<ImputationJob> jobsAfterSave = imputationService.getImputationJobs(testUserId);
    assertEquals(2, jobsAfterSave.size());

    // verify info written to the jobs table
    ImputationJob savedJob = imputationService.getImputationJob(savedUUID, testUserId);
    assertEquals(testJobId, savedJob.getJobId());
    assertEquals(testPipelineVersion, savedJob.getPipelineVersion());
    assertEquals(testUserId, savedJob.getUserId());

    // verify info written to pipelineInputs table
    Optional<PipelineInput> pipelineInput = pipelineInputsRepository.findById(savedJob.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals("{first_key=first_value}", pipelineInput.get().getInputs());
  }

  @Test
  void testWriteDuplicateJob() {
    // try to save a job with the same job id two times, the second time it should throw duplicate
    // exception error
    ImputationJob newJob = createTestJobWithJobId(testJobId);

    ImputationJob savedJobFirst = imputationService.writeJobToDbThrowsDuplicateException(newJob);
    assertNotNull(savedJobFirst);

    ImputationJob newJobSameId = createTestJobWithJobId(testJobId);
    assertThrows(
        DuplicateObjectException.class,
        () -> imputationService.writeJobToDbThrowsDuplicateException(newJobSameId));
  }

  @Test
  void testGetCorrectNumberOfRows() {
    // A test row should exist for this user.
    List<ImputationJob> jobs = imputationJobsRepository.findAllByUserId(testUserId);
    assertEquals(1, jobs.size());

    // insert another row and verify that it shows up
    ImputationJob newJob = createTestJobWithJobId(testJobId);

    imputationJobsRepository.save(newJob);
    jobs = imputationJobsRepository.findAllByUserId(testUserId);
    assertEquals(2, jobs.size());
  }

  @Test
  void testCorrectUserIsolation() {
    // A test row should exist for this user.
    List<ImputationJob> jobs = imputationJobsRepository.findAllByUserId(testUserId);
    assertEquals(1, jobs.size());

    // insert row for second user and verify that it shows up
    String testUserId2 = TestUtils.TEST_USER_ID_2;
    ImputationJob newJob = createTestJobWithJobIdAndUser(UUID.randomUUID(), testUserId2);
    imputationJobsRepository.save(newJob);

    // Verify that the old userid still show only 1 record
    jobs = imputationJobsRepository.findAllByUserId(testUserId);
    assertEquals(1, jobs.size());

    // Verify the new user's id shows a single job as well
    jobs = imputationJobsRepository.findAllByUserId(testUserId2);
    assertEquals(1, jobs.size());
  }
}
