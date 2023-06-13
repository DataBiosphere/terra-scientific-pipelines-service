package bio.terra.pipelines.db;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.DbJob;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.service.JobsService;
import java.time.Instant;
import java.util.List;
import java.util.UUID;
import javax.transaction.Transactional;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DuplicateKeyException;

class JobsDaoTest extends BaseDaoTest {

  @Autowired JobsService jobsService;
  @Autowired JobsRepository jobsRepository;

  private final String testUserId = "testUser";
  private final String testPipelineId = "testPipeline";

  private DbJob createTestJobWithUUID(UUID jobId) {
    return createTestJobWithUUIDAndUser(jobId, testUserId);
  }

  private DbJob createTestJobWithUUIDAndUser(UUID jobId, String userId) {
    Instant timeSubmitted = Instant.now();
    String status = "SUBMITTED";
    return new DbJob(jobId, userId, testPipelineId, "testVersion", timeSubmitted, null, status);
  }

  @Test
  @Transactional
  void testWriteValidJob() {
    String testNewJobId = "deadbeef-dead-beef-aaaa-aaaadeadbeef";
    UUID jobId = UUID.fromString(testNewJobId);
    jobsService.createJob(testUserId, testPipelineId, testNewJobId);

    DbJob savedDbJob = jobsRepository.findByJobId(jobId);
    List<DbJob> dbJobs = (List<DbJob>) jobsRepository.findAll();
    assertEquals(jobId, savedDbJob.getJobId());
  }

  @Test
  void testWriteDuplicateJob() {
    // repeating the exact same logic as above, except using a job id that we know is preloaded in
    // the test data
    String testPresentJobId = "deadbeef-dead-beef-deaf-beefdeadbeef";
    UUID jobId = UUID.fromString(testPresentJobId);
    DbJob newJob = createTestJobWithUUID(jobId);

    // this should not write a job to the db since the job id already exists
    assertThrows(DuplicateKeyException.class, () -> jobsRepository.save(newJob));
  }

  @Test
  void testGetCorrectNumberOfRows() {
    // A test row should exist for this user.
    List<DbJob> jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(1, jobs.size());

    // insert another row and verify that it shows up
    DbJob newJob = createTestJobWithUUID(UUID.randomUUID());

    jobsRepository.save(newJob);
    jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(2, jobs.size());
  }

  @Test
  void testCorrectUserIsolation() {
    // A test row should exist for this user.
    List<DbJob> jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(1, jobs.size());

    // insert row for second user and verify that it shows up
    String testUserId2 = "testUser2";
    DbJob newJob = createTestJobWithUUIDAndUser(UUID.randomUUID(), testUserId2);
    jobsRepository.save(newJob);

    // Verify that the old userid still show only 1 record
    jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId);
    assertEquals(1, jobs.size());

    // Verify the new user's id shows a single job as well
    jobs = jobsRepository.findAllByPipelineIdAndUserId(testPipelineId, testUserId2);
    assertEquals(1, jobs.size());
  }
}
