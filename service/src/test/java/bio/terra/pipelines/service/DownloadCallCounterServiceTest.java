package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.db.entities.DownloadCallCount;
import bio.terra.pipelines.db.repositories.DownloadCallCounterRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class DownloadCallCounterServiceTest extends BaseEmbeddedDbTest {

  @Autowired DownloadCallCounterService downloadCallCounterService;

  @Autowired DownloadCallCounterRepository downloadCallCounterRepository;

  @Test
  void incrementDownloadCallCount() {
    UUID jobId = TestUtils.TEST_NEW_UUID;

    downloadCallCounterService.incrementDownloadCallCount(jobId);

    DownloadCallCount retrievedCallEntity = downloadCallCounterRepository.findByJobId(jobId).get();
    Long id = retrievedCallEntity.getId();
    assertEquals(jobId, retrievedCallEntity.getJobId());
    assertEquals(1, retrievedCallEntity.getDownloadCallCount());
    assertEquals(
        retrievedCallEntity.getFirstCallTimestamp(), retrievedCallEntity.getLatestCallTimestamp());

    // now increment again to verify it updates existing record
    downloadCallCounterService.incrementDownloadCallCount(jobId);

    retrievedCallEntity = downloadCallCounterRepository.findByJobId(jobId).get();
    assertEquals(id, retrievedCallEntity.getId());
    assertEquals(2, retrievedCallEntity.getDownloadCallCount());
    assertTrue(
        retrievedCallEntity
            .getLatestCallTimestamp()
            .isAfter(retrievedCallEntity.getFirstCallTimestamp()));
  }
}
