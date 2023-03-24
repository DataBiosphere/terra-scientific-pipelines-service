package bio.terra.pipelines.db;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.service.model.Job;
import java.time.Instant;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobsDaoTest extends BaseDaoTest {

  @Autowired JobsDao jobsDao;

  private final String testPresentJobId = "deadbeef-dead-beef-deaf-beefdeadbeef";
  private final String testUserId = "testUser";
  private final String testPipelineId = "testPipeline";

  @Test
  void testWriteValidJob() {
    UUID jobId = UUID.randomUUID();
    Instant timeSubmitted = Instant.now();

    String status = "SUBMITTED";
    Job jobFull =
        new Job(jobId, testUserId, testPipelineId, "testVersion", timeSubmitted, null, status);
    UUID uuidResult = jobsDao.createJob(jobFull);

    assertEquals(jobId, uuidResult);
  }

  @Test
  void testWriteDuplicateJob() {
    // repeating the exact same logic as above, except using a job id that we know is preloaded in the test data
    UUID jobId = UUID.fromString(testPresentJobId);
    Instant timeSubmitted = Instant.now();

    String status = "SUBMITTED";
    Job jobFull =
        new Job(jobId, testUserId, testPipelineId, "testVersion", timeSubmitted, null, status);

    //    List<Job> jobs = jobsDao.getJobs(testUserId, testPipelineId);
    assertThrows(DuplicateObjectException.class, () -> jobsDao.createJob(jobFull));
  }

}
