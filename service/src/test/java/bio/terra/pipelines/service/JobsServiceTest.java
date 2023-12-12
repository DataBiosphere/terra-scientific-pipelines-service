package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.time.Instant;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobsServiceTest extends BaseContainerTest {

  @Autowired JobsService jobsService;
  @Autowired JobsRepository jobsRepository;
  @Autowired PipelineInputsRepository pipelineInputsRepository;

  private final String testUserId = "testUser";
  private final String testPipelineId = "testPipeline";
  private final String testPipelineVersion = "testVersion";
  private final Object testPipelineInputs = new LinkedHashMap<>(Map.of("first_key", "first_value"));

  private final String duplicateJobIdString = "deadbeef-dead-beef-aaaa-beefdeadbeef";

  private Job createTestJobWithJobId(UUID jobId) {
    return createTestJobWithJobIdAndUser(jobId, testUserId);
  }

  private Job createTestJobWithJobIdAndUser(UUID jobId, String userId) {
    Instant timeSubmitted = Instant.now();
    String status = "SUBMITTED";
    return new Job(jobId, userId, testPipelineId, testPipelineVersion, timeSubmitted, null, status);
  }

  @Test
  void testWriteJobToDb() {
    List<Job> jobsDefault = jobsService.getJobs(testUserId, testPipelineId);
    // test data migration inserts one row by default
    assertEquals(1, jobsDefault.size());

    String testJobId = jobsService.createJobId();
    String testStatus = "TEST STATUS";
    Instant testTimeSubmitted = Instant.now();

    UUID savedUUID =
        jobsService.writeJobToDb(
            UUID.fromString(testJobId),
            testUserId,
            testPipelineId,
            testPipelineVersion,
            testTimeSubmitted,
            testStatus,
            testPipelineInputs);

    List<Job> jobsAfterSave = jobsService.getJobs(testUserId, testPipelineId);
    assertEquals(2, jobsAfterSave.size());

    // verify info written to the jobs table
    Job savedJob = jobsService.getJob(testUserId, testPipelineId, savedUUID);
    assertEquals(savedUUID, savedJob.getJobId());
    assertEquals(testPipelineId, savedJob.getPipelineId());
    assertEquals(testPipelineVersion, savedJob.getPipelineVersion());
    assertEquals(testUserId, savedJob.getUserId());
    assertEquals(testStatus, savedJob.getStatus());
    assertEquals(testTimeSubmitted, savedJob.getTimeSubmitted());

    // verify info written to pipelineInputs table
    Optional<PipelineInput> pipelineInput = pipelineInputsRepository.findById(savedJob.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals("{first_key=first_value}", pipelineInput.get().getInputs());
  }

  @Test
  void testWriteDuplicateJob() {
    // try to save a job with the same job id two times, the second time it should throw duplicate
    // exception error
    UUID testJobId = UUID.fromString(duplicateJobIdString);

    Job newJob = createTestJobWithJobId(testJobId);

    Job savedJobFirst = jobsService.writeJobToDbThrowsDuplicateException(newJob);
    assertNotNull(savedJobFirst);

    Job newJobSameId = createTestJobWithJobId(testJobId);
    assertThrows(
        DuplicateObjectException.class,
        () -> jobsService.writeJobToDbThrowsDuplicateException(newJobSameId));
  }

  @Test
  void testGetCorrectNumberOfRows() {
    // A test row should exist for this user.
    List<Job> jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(1, jobs.size());

    // insert another row and verify that it shows up
    Job newJob = createTestJobWithJobId(UUID.randomUUID());

    jobsRepository.save(newJob);
    jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(2, jobs.size());
  }

  @Test
  void testCorrectUserIsolation() {
    // A test row should exist for this user.
    List<Job> jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(1, jobs.size());

    // insert row for second user and verify that it shows up
    String testUserId2 = "testUser2";
    Job newJob = createTestJobWithJobIdAndUser(UUID.randomUUID(), testUserId2);
    jobsRepository.save(newJob);

    // Verify that the old userid still show only 1 record
    jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(1, jobs.size());

    // Verify the new user's id shows a single job as well
    jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId2);
    assertEquals(1, jobs.size());
  }

  @Test
  void testGetTimestamp() {
    Instant timestamp1 = jobsService.getCurrentTimestamp();
    Instant timestamp2 = jobsService.getCurrentTimestamp();

    assertNotNull(timestamp1);
    assertNotNull(timestamp2);
    assertNotEquals(timestamp1, timestamp2);
    // timestamp1 should be before timestamp2
    assertTrue(timestamp1.compareTo(timestamp2) < 0);
  }

  @Test
  void testCreateJobId() {
    String jobId1 = jobsService.createJobId();
    String jobId2 = jobsService.createJobId();

    assertNotNull(jobId1);
    assertNotNull(jobId2);
    assertNotEquals(jobId1, jobId2);

    // make sure the job id is a valid UUID
    assertEquals(jobId1, UUID.fromString(jobId1).toString());
    assertEquals(jobId2, UUID.fromString(jobId2).toString());
  }
}
