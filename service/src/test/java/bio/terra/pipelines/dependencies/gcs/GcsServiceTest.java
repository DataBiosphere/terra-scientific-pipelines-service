package bio.terra.pipelines.dependencies.gcs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.times;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.cloud.storage.*;
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
  @Captor private ArgumentCaptor<Storage.CopyRequest> copyRequestCaptor;
  @Captor private ArgumentCaptor<BlobId> blobIdCaptor;
  private final String bucketName = "bucketName";
  private final String objectName = "objectName";
  private final GcsFile gcsFile = new GcsFile("gs://%s/%s".formatted(bucketName, objectName));
  private final String userBearerToken = "userBearerToken";

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

    when(gcsClient.getStorageService()).thenReturn(mockStorageService);
    when(gcsClient.getStorageService(null)).thenReturn(mockStorageService);
    when(gcsClient.getStorageService(userBearerToken)).thenReturn(mockStorageService);
  }

  @Test
  void serviceHasBlobReadAccessTrue() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);
    Blob mockBlob = mock(Blob.class);
    when(mockBlob.exists()).thenReturn(true);

    when(mockStorageService.get(blobId, blobOption)).thenReturn(mockBlob);
    assertTrue(gcsService.serviceHasFileReadAccess(gcsFile));
  }

  @Test
  void serviceHasBlobReadAccessFalse() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);
    Blob mockBlob = mock(Blob.class);
    when(mockBlob.exists()).thenReturn(false);

    when(mockStorageService.get(blobId, blobOption)).thenReturn(mockBlob);
    assertFalse(gcsService.serviceHasFileReadAccess(gcsFile));
  }

  @Test
  void serviceHasBlobReadAccessFalseNoBlob() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);

    when(mockStorageService.get(blobId, blobOption)).thenReturn(null);
    assertFalse(gcsService.serviceHasFileReadAccess(gcsFile));
  }

  @Test
  void serviceHasBlobReadAccessFalseStorageException() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);

    when(mockStorageService.get(blobId, blobOption))
        .thenThrow(new StorageException(500, "Storage exception"));
    assertFalse(gcsService.serviceHasFileReadAccess(gcsFile));
  }

  @Test
  void userHasBlobReadAccessTrue() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);
    Blob mockBlob = mock(Blob.class);
    when(mockBlob.exists()).thenReturn(true);

    when(mockStorageService.get(blobId, blobOption)).thenReturn(mockBlob);
    assertTrue(gcsService.userHasFileReadAccess(gcsFile, userBearerToken));
  }

  @Test
  void userHasBlobReadAccessFalse() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);
    Blob mockBlob = mock(Blob.class);
    when(mockBlob.exists()).thenReturn(false);

    when(mockStorageService.get(blobId, blobOption)).thenReturn(mockBlob);
    assertFalse(gcsService.userHasFileReadAccess(gcsFile, userBearerToken));
  }

  @Test
  void userHasBlobReadAccessFalseNoBlob() {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    Storage.BlobGetOption blobOption = Storage.BlobGetOption.fields(Storage.BlobField.NAME);

    when(mockStorageService.get(blobId, blobOption)).thenReturn(null);
    assertFalse(gcsService.userHasFileReadAccess(gcsFile, userBearerToken));
  }

  @Test
  void userHasBlobReadAccessFalseNullToken() {
    assertThrows(NullPointerException.class, () -> gcsService.userHasFileReadAccess(gcsFile, null));
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
            any(Storage.SignUrlOption.class)))
        .thenReturn(fakeURL);

    URL generatedURL = gcsService.generatePutObjectSignedUrl(bucketName, objectName);

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

    URL generatedURL = gcsService.generateResumablePostObjectSignedUrl(bucketName, objectName);

    assertEquals(fakeURL, generatedURL);

    BlobInfo blobInfo = blobInfoCaptor.getValue();
    assertEquals(bucketName, blobInfo.getBucket());
    assertEquals(objectName, blobInfo.getName());
  }

  @Test
  void generateGetObjectSignedUrl() {
    URL fakeSignedUrl = getFakeURL();
    when(mockStorageService.signUrl(
            blobInfoCaptor.capture(),
            eq(testSignedUrlGetDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenReturn(fakeSignedUrl);

    URL generatedSignedUrl = gcsService.generateGetObjectSignedUrl(gcsFile);
    assertEquals(fakeSignedUrl, generatedSignedUrl);

    BlobInfo blobInfo = blobInfoCaptor.getValue();
    assertEquals(bucketName, blobInfo.getBucket());
    assertEquals(objectName, blobInfo.getName());
  }

  @Test
  void generateGetObjectSignedUrlFileWithSlashes() {
    URL fakeURL = getFakeURL();
    when(mockStorageService.signUrl(
            blobInfoCaptor.capture(),
            eq(testSignedUrlGetDuration),
            eq(TimeUnit.HOURS),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class),
            any(Storage.SignUrlOption.class)))
        .thenReturn(fakeURL);

    String objectNameWithSlash = "objectName/with/slash";
    GcsFile gcsFileForSignedUrl =
        new GcsFile("gs://%s/%s".formatted(bucketName, objectNameWithSlash));
    URL generatedURL = gcsService.generateGetObjectSignedUrl(gcsFileForSignedUrl);
    assertEquals(fakeURL, generatedURL);

    BlobInfo blobInfo = blobInfoCaptor.getValue();
    assertEquals(bucketName, blobInfo.getBucket());
    assertEquals(objectNameWithSlash, blobInfo.getName());
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

    URL generatedURL = gcsService.generatePutObjectSignedUrl(bucketName, objectName);
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
          gcsService.generatePutObjectSignedUrl(bucketName, objectName);
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
          gcsService.generatePutObjectSignedUrl(bucketName, objectName);
        });
  }

  @Test
  void copyObjectSuccess() {
    GcsFile sourceGcsPath = new GcsFile("gs://source-bucket/path/to/file.vcf.gz");
    GcsFile destinationGcsPath = new GcsFile("gs://destination-bucket/jobId/file.vcf.gz");

    Blob mockBlob = mock(Blob.class);
    CopyWriter mockCopyWriter = mock(CopyWriter.class);

    when(mockCopyWriter.getResult()).thenReturn(mockBlob);

    when(mockStorageService.copy(any(Storage.CopyRequest.class))).thenReturn(mockCopyWriter);
    when(mockStorageService.get(any(BlobId.class))).thenReturn(mockBlob);

    gcsService.copyObject(sourceGcsPath, destinationGcsPath);

    // Verify copy and get were called
    verify(mockStorageService, times(1)).copy(any(Storage.CopyRequest.class));
    verify(mockStorageService, times(1)).get(any(BlobId.class));
  }

  @Test
  void copyObjectRetriesEventuallySucceed() {
    GcsFile sourceGcsPath = new GcsFile("gs://source-bucket/path/to/file.vcf.gz");
    GcsFile destinationPath = new GcsFile("gs://destination-bucket/jobId/file.vcf.gz");

    Blob mockBlob = mock(Blob.class);
    CopyWriter mockCopyWriter = mock(CopyWriter.class);

    when(mockCopyWriter.getResult()).thenReturn(mockBlob);

    when(mockStorageService.copy(any(Storage.CopyRequest.class)))
        .thenAnswer(errorAnswer) // first call fails
        .thenReturn(mockCopyWriter); // retry succeeds

    when(mockStorageService.get(any(BlobId.class))).thenReturn(mockBlob);

    gcsService.copyObject(sourceGcsPath, destinationPath);

    // Verify copy was called twice (once failed, once succeeded) and get was called once
    verify(mockStorageService, times(2)).copy(any(Storage.CopyRequest.class));
    verify(mockStorageService, times(1)).get(any(BlobId.class));
  }

  @Test
  void copyObjectStorageExceptionDoNotRetry() {
    GcsFile sourceGcsPath = new GcsFile("gs://source-bucket/path/to/file.vcf.gz");
    GcsFile destinationPath = new GcsFile("gs://destination-bucket/jobId/file.vcf.gz");

    when(mockStorageService.copy(any(Storage.CopyRequest.class)))
        .thenThrow(new StorageException(400, "Storage exception"));

    assertThrows(
        GcsServiceException.class,
        () -> {
          gcsService.copyObject(sourceGcsPath, destinationPath);
        });
  }

  @Test
  void deleteObjectSuccess() {
    String blobName = "path/to/file.vcf.gz";

    when(mockStorageService.delete(any(BlobId.class))).thenReturn(true);

    gcsService.deleteObject(bucketName, blobName);

    // Verify delete was called once
    verify(mockStorageService, times(1)).delete(any(BlobId.class));
  }

  @Test
  void deleteObjectNotFound() {
    String blobName = "path/to/nonexistent-file.vcf.gz";

    when(mockStorageService.delete(any(BlobId.class))).thenReturn(false);

    gcsService.deleteObject(bucketName, blobName);

    // Verify delete was called once even though object was not found
    verify(mockStorageService, times(1)).delete(any(BlobId.class));
  }

  @Test
  void deleteObjectSocketExceptionRetriesEventuallySucceed() {
    String blobName = "path/to/file.vcf.gz";

    when(mockStorageService.delete(any(com.google.cloud.storage.BlobId.class)))
        .thenAnswer(errorAnswer)
        .thenReturn(true);

    gcsService.deleteObject(bucketName, blobName);

    // Verify delete was called twice (once failed, once succeeded)
    verify(mockStorageService, times(2)).delete(any(BlobId.class));
  }

  @Test
  void deleteObjectSocketExceptionRetriesEventuallyFail() {
    String blobName = "path/to/file.vcf.gz";

    when(mockStorageService.delete(any(BlobId.class)))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer);

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          gcsService.deleteObject(bucketName, blobName);
        });
  }

  @Test
  void deleteObjectStorageExceptionDoNotRetry() {
    String blobName = "path/to/file.vcf.gz";

    when(mockStorageService.delete(any(BlobId.class)))
        .thenThrow(new StorageException(400, "Storage exception"));

    assertThrows(
        GcsServiceException.class,
        () -> {
          gcsService.deleteObject(bucketName, blobName);
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
