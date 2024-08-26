package bio.terra.pipelines.dependencies.gcs;

import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.cloud.storage.Storage;
import java.net.MalformedURLException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class GcsServiceTest extends BaseEmbeddedDbTest {

  @Autowired @InjectMocks private GcsService gcsService;
  @MockBean private GcsClient gcsClient;

  private Storage mockStorageService = mock(Storage.class);

  @BeforeEach
  void setup() {
    when(gcsClient.getStorageService(any())).thenReturn(mockStorageService);
  }

  @Test
  void generatePutObjectSignedUrl() throws MalformedURLException {
    //    URL fakeURL = new URL("https://storage.googleapis.com/signed-url-stuff");
    //    when(mockStorageService.signUrl(
    //            any(BlobInfo.class), anyLong(), any(TimeUnit.class),
    // any(Storage.SignUrlOption.class)))
    //        .thenReturn(fakeURL);
    //
    //    URL generatedURL =
    //        gcsService.generatePutObjectSignedUrl("projectId", "bucketName", "objectName");
    //    assertEquals(fakeURL, generatedURL);
  }
}
