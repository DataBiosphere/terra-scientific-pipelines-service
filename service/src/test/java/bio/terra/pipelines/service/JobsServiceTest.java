package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.time.Instant;
import java.util.List;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobsServiceTest extends BaseContainerTest {

  @Autowired JobsService jobsService;
  @Autowired JobsRepository jobsRepository;

  private final String testUserId = "testUser";
  private final String testPipelineId = "testPipeline";
  private final String testPipelineVersion = "testVersion";

  private Job createTestJobWithJobId(UUID jobId) {
    return createTestJobWithJobIdAndUser(jobId, testUserId);
  }

  private Job createTestJobWithJobIdAndUser(UUID jobId, String userId) {
    Instant timeSubmitted = Instant.now();
    String status = "SUBMITTED";
    return new Job(jobId, userId, testPipelineId, testPipelineVersion, timeSubmitted, null, status);
  }

  @Test
  void testWriteValidJob() {
    List<Job> jobsDefault = jobsService.getJobs(testUserId, testPipelineId);
    // test data migration inserts one row by default
    assertEquals(1, jobsDefault.size());

    UUID savedUUID = jobsService.createJob(testUserId, testPipelineId, testPipelineVersion);

    List<Job> jobsAfterSave = jobsService.getJobs(testUserId, testPipelineId);
    assertEquals(2, jobsAfterSave.size());

    Job savedJob = jobsService.getJob(testUserId, testPipelineId, savedUUID);
    assertEquals(savedJob.getJobId(), savedUUID);
    assertEquals(savedJob.getPipelineId(), testPipelineId);
    assertEquals(savedJob.getPipelineVersion(), testPipelineVersion);
    assertEquals(savedJob.getUserId(), testUserId);
  }

  @Test
  void testWriteDuplicateJob() {
    // try to save a job with the same job id two times, the second time it should not save and
    // return null
    UUID testJobId = UUID.fromString("deadbeef-dead-beef-aaaa-beefdeadbeef");

    Job newJob = createTestJobWithJobId(testJobId);

    UUID savedJobUUIDFirst = jobsService.writeJobToDbRetryDuplicateException(newJob);
    assertNotNull(savedJobUUIDFirst);
    Job jobWithDuplicateJobId = createTestJobWithJobId(testJobId);
    UUID savedJobUUIDSecond =
        jobsService.writeJobToDbRetryDuplicateException(jobWithDuplicateJobId);
    // this should not write a job to the db since the job id already exists and thus will return
    // null
    assertNull(savedJobUUIDSecond);
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
}
