package bio.terra.pipelines.db;

import bio.terra.pipelines.service.model.Job;
import java.time.Instant;
import java.util.UUID;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class JobsDaoTest extends BaseDaoTest {

  @Autowired JobsDao jobsDao;

  private final String testUserId = "testUser";
  private final String testPipelineId = "testPipeline";

  @Test
  void testWriteJob() {
    UUID jobId = UUID.randomUUID();
    Instant timeSubmitted = Instant.now();

    String status = "SUBMITTED";
    Job jobFull =
        new Job(jobId, testUserId, testPipelineId, "testVersion", timeSubmitted, null, status);
    UUID uuidResult = jobsDao.createJob(jobFull);

    assertEquals(jobId, uuidResult);
  }
}
