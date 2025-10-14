package bio.terra.pipelines.dependencies.gcs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.cloud.storage.BlobInfo;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageException;
import java.net.MalformedURLException;
import java.net.SocketTimeoutException;
import java.net.URL;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Captor;
import org.mockito.InjectMocks;
import org.mockito.stubbing.Answer;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class GcsServiceTest extends BaseEmbeddedDbTest {

  @Autowired @InjectMocks private GcsService gcsService;
  @MockitoBean private GcsClient gcsClient;

  private final Storage mockStorageService = mock(Storage.class);

  @Captor private ArgumentCaptor<BlobInfo> blobInfoCaptor;
  private final String projectId = "projectId";
  private final String bucketName = "bucketName";
  private final String objectName = "objectName";

  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  private final Long testSignedUrlPutDuration = 8L;
  private final Long testSignedUrlGetDuration = 1L;

  final Answer<Object> errorAnswer =
      invocation -> {
        throw new SocketTimeoutException("Timeout");
      };

  @BeforeEach
  void setup() {
    FixedBackOffPolicy smallerBackoff = new FixedBackOffPolicy();
    smallerBackoff.setBackOffPeriod(5L); // 5 ms
    template.setBackOffPolicy(smallerBackoff);

    when(gcsClient.getStorageService(projectId)).thenReturn(mockStorageService);
  }

  private URL getFakeURL() {
    try {
      return new URL("https://storage.googleapis.com/signed-url-stuff?X-Goog-Signature=12345");
    } catch (MalformedURLException e) {
      return null;
    }
  }

  @Test
  void generatePutObjectSignedUrl() {
    URL fakeURL = getFakeURL();
    when(mockStorageService.signUrl(
            blobInfoCaptor.capture(),
            eq(testSignedUrlPutDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenReturn(fakeURL);

    URL generatedURL = gcsService.generatePutObjectSignedUrl(projectId, bucketName, objectName);

    assertEquals(fakeURL, generatedURL);

    BlobInfo blobInfo = blobInfoCaptor.getValue();
    assertEquals(bucketName, blobInfo.getBucket());
    assertEquals(objectName, blobInfo.getName());
  }

  @Test
  void generateResumablePostObjectSignedUrl() {
    URL fakeURL = getFakeURL();
    when(mockStorageService.signUrl(
            blobInfoCaptor.capture(),
            eq(testSignedUrlPutDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenReturn(fakeURL);

    URL generatedURL =
        gcsService.generateResumablePostObjectSignedUrl(projectId, bucketName, objectName);

    assertEquals(fakeURL, generatedURL);

    BlobInfo blobInfo = blobInfoCaptor.getValue();
    assertEquals(bucketName, blobInfo.getBucket());
    assertEquals(objectName, blobInfo.getName());
  }

  @Test
  void generateGetObjectSignedUrl() {
    URL fakeURL = getFakeURL();
    when(mockStorageService.signUrl(
            blobInfoCaptor.capture(),
            eq(testSignedUrlGetDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenReturn(fakeURL);

    URL generatedURL = gcsService.generateGetObjectSignedUrl(projectId, bucketName, objectName);
    assertEquals(fakeURL, generatedURL);

    BlobInfo blobInfo = blobInfoCaptor.getValue();
    assertEquals(bucketName, blobInfo.getBucket());
    assertEquals(objectName, blobInfo.getName());
  }

  @Test
  void socketExceptionRetriesEventuallySucceed() {
    URL fakeURL = getFakeURL();

    when(mockStorageService.signUrl(
            any(BlobInfo.class),
            eq(testSignedUrlPutDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenAnswer(errorAnswer)
        .thenReturn(fakeURL);

    URL generatedURL = gcsService.generatePutObjectSignedUrl(projectId, bucketName, objectName);
    assertEquals(fakeURL, generatedURL);
  }

  @Test
  void socketExceptionRetriesEventuallyFail() {
    URL fakeURL = getFakeURL();

    when(mockStorageService.signUrl(
            any(BlobInfo.class),
            eq(testSignedUrlPutDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(fakeURL);

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          gcsService.generatePutObjectSignedUrl(projectId, bucketName, objectName);
        });
  }

  @Test
  void storageExceptionDoNotRetry() {
    when(mockStorageService.signUrl(
            any(BlobInfo.class),
            eq(testSignedUrlPutDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenThrow(new StorageException(400, "Storage exception"));

    assertThrows(
        GcsServiceException.class,
        () -> {
          gcsService.generatePutObjectSignedUrl(projectId, bucketName, objectName);
        });
  }

  @Test
  void cleanSignedUrl() throws MalformedURLException {
    // signed URL with X-Goog-Signature as last element
    URL fakeURLSignatureLast =
        new URL(
            "https://storage.googleapis.com/fc-secure-6970c3a9-dc92-436d-af3d-917bcb4cf05a/user-input-files/ffaffa12-5717-4562-b3fc-2c963f66afa6/TEST.vcf.gz?X-Goog-Date=20240823T170006Z&X-Goog-Expires=900&X-Goog-SignedHeaders=content-type%3Bhost&X-Goog-Signature=12345");
    String expectedCleanedURLSignatureLast =
        "https://storage.googleapis.com/fc-secure-6970c3a9-dc92-436d-af3d-917bcb4cf05a/user-input-files/ffaffa12-5717-4562-b3fc-2c963f66afa6/TEST.vcf.gz?X-Goog-Date=20240823T170006Z&X-Goog-Expires=900&X-Goog-SignedHeaders=content-type%3Bhost&X-Goog-Signature=REDACTED";
    assertEquals(expectedCleanedURLSignatureLast, GcsService.cleanSignedUrl(fakeURLSignatureLast));

    // signed URL with no X-Goog-Signature should not throw an exception
    URL fakeURLNoSignature =
        new URL(
            "https://storage.googleapis.com/fc-secure-6970c3a9-dc92-436d-af3d-917bcb4cf05a/user-input-files/ffaffa12-5717-4562-b3fc-2c963f66afa6/TEST.vcf.gz?X-Goog-Date=20240823T170006Z&X-Goog-Expires=900");
    assertEquals(fakeURLNoSignature.toString(), GcsService.cleanSignedUrl(fakeURLNoSignature));

    // signed URL with X-Goog-Signature in the middle
    URL fakeURLSignatureMiddle =
        new URL(
            "https://storage.googleapis.com/fc-secure-6970c3a9-dc92-436d-af3d-917bcb4cf05a/user-input-files/ffaffa12-5717-4562-b3fc-2c963f66afa6/TEST.vcf.gz?X-Goog-Date=20240823T170006Z&X-Goog-Expires=900&X-Goog-Signature=12345&X-Goog-SignedHeaders=content-type%3Bhost&Last-Element=foobar");
    String expectedCleanedURLSignatureMiddle =
        "https://storage.googleapis.com/fc-secure-6970c3a9-dc92-436d-af3d-917bcb4cf05a/user-input-files/ffaffa12-5717-4562-b3fc-2c963f66afa6/TEST.vcf.gz?X-Goog-Date=20240823T170006Z&X-Goog-Expires=900&X-Goog-SignedHeaders=content-type%3Bhost&Last-Element=foobar&X-Goog-Signature=REDACTED";
    assertEquals(
        expectedCleanedURLSignatureMiddle, GcsService.cleanSignedUrl(fakeURLSignatureMiddle));
  }
}
