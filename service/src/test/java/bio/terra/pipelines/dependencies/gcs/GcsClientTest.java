package bio.terra.pipelines.dependencies.gcs;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.cloud.storage.Storage;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class GcsClientTest extends BaseEmbeddedDbTest {
  @Autowired GcsClient gcsClient;

  @Test
  void getGcsStorageService() {
    String projectId = "test-project-id";
    Storage storageService = gcsClient.getStorageService(projectId);

    assertEquals(projectId, storageService.getOptions().getProjectId());
  }
}
