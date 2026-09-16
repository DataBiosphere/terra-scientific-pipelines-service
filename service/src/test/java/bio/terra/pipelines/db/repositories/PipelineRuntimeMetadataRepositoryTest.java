package bio.terra.pipelines.db.repositories;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.time.Instant;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelineRuntimeMetadataRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;

  @Test
  // This test verifies that the pipeline_runtime_metadata_updated_trigger correctly updates the
  // `updated` column to the current time whenever a row is modified.
  //
  // Note: this test is intentionally NOT @Transactional. Postgres' now() returns the start time
  // of the enclosing transaction, so the insert and update below must run in separate
  // transactions in order to observe the `updated` column actually change. Each repository call
  // below runs in its own transaction. Test data is cleaned up by the embedded DB refresh that
  // runs after each test method.
  void updatingRowUpdatesUpdatedTimestamp() {
    String pipelineKey = buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, 99);
    PipelineRuntimeMetadata meta = new PipelineRuntimeMetadata(pipelineKey);
    meta.setHidden(true);
    pipelineRuntimeMetadataRepository.save(meta);

    PipelineRuntimeMetadata created =
        pipelineRuntimeMetadataRepository.findById(pipelineKey).orElseThrow();
    Instant createdTimestamp = created.getUpdated();
    assertNotNull(createdTimestamp);

    created.setToolVersion("1.2.3");
    pipelineRuntimeMetadataRepository.save(created);

    PipelineRuntimeMetadata updated =
        pipelineRuntimeMetadataRepository.findById(pipelineKey).orElseThrow();
    assertEquals("1.2.3", updated.getToolVersion());
    assertNotNull(updated.getUpdated());
    assertTrue(
        updated.getUpdated().isAfter(createdTimestamp),
        "expected updated timestamp %s to be after original timestamp %s"
            .formatted(updated.getUpdated(), createdTimestamp));
  }
}
