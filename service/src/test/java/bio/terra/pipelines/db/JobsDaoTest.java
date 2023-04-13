package bio.terra.pipelines.db;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.service.model.Job;
import java.time.Instant;
import java.util.List;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobsDaoTest extends BaseDaoTest {

  @Autowired JobsDao jobsDao;

  private final String testUserId = "testUser";
  private final String testPipelineId = "testPipeline";

  private Job createTestJobWithUUID(UUID jobId) {
    return createTestJobWithUUIDAndUser(jobId, testUserId);
  }

  private Job createTestJobWithUUIDAndUser(UUID jobId, String userId) {
    Instant timeSubmitted = Instant.now();
    String status = "SUBMITTED";
    return new Job(jobId, userId, testPipelineId, "testVersion", timeSubmitted, null, status);
  }

  @Test
  void testWriteValidJob() {
    String testNewJobId = "deadbeef-dead-beef-aaaa-aaaadeadbeef";
    UUID jobId = UUID.fromString(testNewJobId);
    Job newJob = createTestJobWithUUID(jobId);
    UUID uuidResult = jobsDao.createJob(newJob);

    assertEquals(jobId, uuidResult);
  }

  @Test
  void testWriteDuplicateJob() {
    // repeating the exact same logic as above, except using a job id that we know is preloaded in
    // the test data
    String testPresentJobId = "deadbeef-dead-beef-deaf-beefdeadbeef";
    UUID jobId = UUID.fromString(testPresentJobId);
    Job newJob = createTestJobWithUUID(jobId);

    // this should not write a job to the db since the job id already exists
    UUID uuidResult = jobsDao.createJob(newJob);

    assertNull(uuidResult);
  }

  @Test
  void testGetCorrectNumberOfRows() {
    // A test row should exist for this user.
    List<Job> jobs = jobsDao.getJobs(testUserId, testPipelineId);
    assertEquals(1, jobs.size());

    // insert another row and verify that it shows up
    Job newJob = createTestJobWithUUID(UUID.randomUUID());

    jobsDao.createJob(newJob);
    jobs = jobsDao.getJobs(testUserId, testPipelineId);
    assertEquals(2, jobs.size());
  }

  @Test
  void testCorrectUserIsolation() {
    // A test row should exist for this user.
    List<Job> jobs = jobsDao.getJobs(testUserId, testPipelineId);
    assertEquals(1, jobs.size());

    // insert row for second user and verify that it shows up
    String testUserId2 = "testUser2";
    Job newJob = createTestJobWithUUIDAndUser(UUID.randomUUID(), testUserId2);
    jobsDao.createJob(newJob);

    // Verify that the old userid still show only 1 record
    jobs = jobsDao.getJobs(testUserId, testPipelineId);
    assertEquals(1, jobs.size());

    // Verify the new user's id shows a single job as well
    jobs = jobsDao.getJobs(testUserId2, testPipelineId);
    assertEquals(1, jobs.size());
  }
}
