package bio.terra.pipelines.dependencies.gcs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.cloud.storage.Storage;
import java.io.IOException;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class GcsClientTest extends BaseEmbeddedDbTest {
  @Autowired GcsClient gcsClient;

  @Test
  void getGcsStorageServiceUserCreds() throws IOException {
    String userToken = "user-token";
    Storage storageService = gcsClient.getStorageService(userToken);

    assertEquals(
        "Bearer %s".formatted(userToken),
        storageService
            .getOptions()
            .getCredentials()
            .getRequestMetadata()
            .get("Authorization")
            .get(0));
  }

  @Test
  void getGcsStorageServiceDefaultCreds() {
    Storage storageService = gcsClient.getStorageService();

    // The credentials should be null, since default credentials will be added by the Google Cloud
    // client library when making requests
    assertNull(storageService.getOptions().getCredentials());
  }

  @Test
  void getGcsStorageServiceDefaultCredsNull() {
    Storage storageService = gcsClient.getStorageService(null);

    // The credentials should be null, since default credentials will be added by the Google Cloud
    // client library when making requests
    assertNull(storageService.getOptions().getCredentials());
  }

  @Test
  void getGcsStorageServiceDefaultCredsEmpty() {
    Storage storageService = gcsClient.getStorageService("");

    // The credentials should be null, since default credentials will be added by the Google Cloud
    // client library when making requests
    assertNull(storageService.getOptions().getCredentials());
  }
}
