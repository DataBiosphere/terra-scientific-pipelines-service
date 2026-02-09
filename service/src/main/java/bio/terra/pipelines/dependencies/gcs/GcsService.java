package bio.terra.pipelines.dependencies.gcs;

import bio.terra.pipelines.app.configuration.external.GcsConfiguration;
import com.google.cloud.storage.Blob;
import com.google.cloud.storage.BlobId;
import com.google.cloud.storage.BlobInfo;
import com.google.cloud.storage.HttpMethod;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageException;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

/** class to encapsulate interacting with GCS client */
@Service
public class GcsService {

  private final GcsClient gcsClient;
  private final GcsConfiguration gcsConfiguration;
  private final RetryTemplate listenerResetRetryTemplate;

  private static final Logger logger = LoggerFactory.getLogger(GcsService.class);

  public GcsService(
      GcsClient gcsClient,
      GcsConfiguration gcsConfiguration,
      RetryTemplate listenerResetRetryTemplate) {
    this.gcsClient = gcsClient;
    this.gcsConfiguration = gcsConfiguration;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  private final String READ_PERMISSION = "storage.objects.get";
  private final String WRITE_PERMISSION = "storage.objects.create";

  public boolean serviceHasBucketReadAccess(String bucketName) {
    boolean hasAccess = hasBucketPermission(bucketName, READ_PERMISSION, null);
    logAccessCheckResult(
        "Teaspoons service account", READ_PERMISSION, "bucket %s".formatted(bucketName), hasAccess);
    return hasAccess;
  }

  public boolean userHasBucketReadAccess(String bucketName, String accessToken) {
    if (accessToken == null) {
      logger.error("Access token is required to check user bucket read access");
      return false;
    }
    boolean hasAccess = hasBucketPermission(bucketName, READ_PERMISSION, accessToken);
    logAccessCheckResult("User", READ_PERMISSION, "bucket %s".formatted(bucketName), hasAccess);
    return hasAccess;
  }

  public boolean serviceHasBucketWriteAccess(String bucketName) {
    boolean hasAccess = hasBucketPermission(bucketName, WRITE_PERMISSION, null);
    logAccessCheckResult(
        "Teaspoons service account",
        WRITE_PERMISSION,
        "bucket %s".formatted(bucketName),
        hasAccess);
    return hasAccess;
  }

  public boolean userHasBucketWriteAccess(String bucketName, String accessToken) {
    if (accessToken == null) {
      logger.error("Access token is required to check user bucket write access");
      return false;
    }
    boolean hasAccess = hasBucketPermission(bucketName, WRITE_PERMISSION, accessToken);
    logAccessCheckResult("User", WRITE_PERMISSION, "bucket %s".formatted(bucketName), hasAccess);
    return hasAccess;
  }

  /**
   * Helper method to log the result of an access check in a consistent format. Logs at INFO level
   * for granted access and ERROR level for denied access.
   */
  private void logAccessCheckResult(
      String subject, String permission, String resource, boolean hasAccess) {
    if (hasAccess) {
      logger.info("{} has {} access on {}", subject, permission, resource);
    } else {
      logger.error("{} does not have {} access on {}", subject, permission, resource);
    }
  }

  /**
   * Check if the service account has a particular access permission on the bucket.
   *
   * @param bucketName without a prefix
   * @return true if the permission is granted
   */
  private boolean hasBucketPermission(String bucketName, String permission, String accessToken)
      throws StorageException {

    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> {
          List<Boolean> accessResult =
              gcsClient
                  .getStorageService(accessToken)
                  .testIamPermissions(bucketName, List.of(permission));
          return accessResult.size() == 1 ? accessResult.get(0) : false;
        });
  }

  public boolean serviceHasBlobReadAccess(String blobPath) {
    boolean hasAccess = hasBlobReadAccess(blobPath, null);
    logAccessCheckResult(
        "Teaspoons service account", READ_PERMISSION, "blob %s".formatted(blobPath), hasAccess);
    return hasAccess;
  }

  public boolean userHasBlobReadAccess(String blobPath, String accessToken) {
    if (accessToken == null) {
      logger.error("Access token is required to check user blob read access");
      return false;
    }
    boolean hasAccess = hasBlobReadAccess(blobPath, accessToken);
    logAccessCheckResult("User", READ_PERMISSION, "blob %s".formatted(blobPath), hasAccess);
    return hasAccess;
  }

  /**
   * Check if a given user has read access to a GCS file.
   *
   * @param blobPath fully qualified GCS path to the file, e.g. `gs://bucket/path/to/file.txt`
   * @param accessToken for the calling user. pass null to use application default credentials.
   * @return boolean whether the caller has read access to the GCS file
   */
  private boolean hasBlobReadAccess(String blobPath, String accessToken) {
    BlobId blobId = BlobId.fromGsUtilUri(blobPath);
    try {
      // Attempt to retrieve a client-side representation of the blob with minimal metadata
      Blob blob =
          gcsClient
              .getStorageService(accessToken)
              .get(blobId, Storage.BlobGetOption.fields(Storage.BlobField.NAME));

      if (blob != null && blob.exists()) {
        return true;
      }
    } catch (StorageException e) {
      logger.error("An error occurred checking blob access: " + e.getMessage());
      return false;
    }
    return false;
  }

  /**
   * Generates and returns a PUT (write-only) signed url for a specific object in a bucket. See
   * documentation on signed urls <a
   * href="https://cloud.google.com/storage/docs/access-control/signed-urls">here</a>.
   *
   * <p>The output URL can be used with a curl command to upload an object to the destination: `curl
   * -X PUT -H 'Content-Type: application/octet-stream' --upload-file my-file '{url}'`
   *
   * @param projectId Google project id
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object (subdirectories + file name)
   * @return url that can be used to write an object to GCS
   */
  public URL generatePutObjectSignedUrl(String projectId, String bucketName, String objectName)
      throws StorageException {
    return generateUploadSignedUrl(projectId, bucketName, objectName, false);
  }

  /**
   * Generates and returns a resumable POST (write-only) signed url for a specific object in a
   * bucket. See documentation on signed urls <a
   * href="https://cloud.google.com/storage/docs/access-control/signed-urls">here</a> and resumable
   * uploads <a
   * href="https://cloud.google.com/storage/docs/performing-resumable-uploads#initiate-session">here</a>.
   *
   * @param projectId Google project id
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object (subdirectories + file name)
   * @return url that can be used to initiate a resumable object upload to GCS
   */
  public URL generateResumablePostObjectSignedUrl(
      String projectId, String bucketName, String objectName) throws StorageException {
    return generateUploadSignedUrl(projectId, bucketName, objectName, true);
  }

  private URL generateUploadSignedUrl(
      String projectId, String bucketName, String objectName, boolean isResumable)
      throws StorageException {
    // define target blob object resource
    BlobInfo blobInfo = BlobInfo.newBuilder(BlobId.of(bucketName, objectName)).build();

    // generate signed URL
    Map<String, String> extensionHeaders = new HashMap<>();
    extensionHeaders.put("Content-Type", "application/octet-stream");
    if (isResumable) {
      extensionHeaders.put("x-goog-resumable", "start");
    }

    URL url =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                gcsClient
                    .getStorageServiceWithProject(projectId)
                    .signUrl(
                        blobInfo,
                        gcsConfiguration.signedUrlPutDurationHours(),
                        TimeUnit.HOURS,
                        Storage.SignUrlOption.httpMethod(
                            isResumable ? HttpMethod.POST : HttpMethod.PUT),
                        Storage.SignUrlOption.withExtHeaders(extensionHeaders),
                        Storage.SignUrlOption.withV4Signature()));
    String cleanSignedUrlString = cleanSignedUrl(url);
    logger.info(
        "Generated {} signed URL: {}",
        isResumable ? "resumable POST" : "PUT",
        cleanSignedUrlString);

    return url;
  }

  /**
   * Generates and returns a GET (read-only) signed url for a specific object in a bucket. See
   * documentation on signed urls <a
   * href="https://cloud.google.com/storage/docs/access-control/signed-urls">here</a>.
   *
   * <p>The output URL can be used with a curl command to download an object: `curl '{url}' >
   * {local_file_name}`
   *
   * @param projectId Google project id
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object (subdirectories + file name)
   * @return url that can be used to download an object to GCS
   */
  public URL generateGetObjectSignedUrl(String projectId, String bucketName, String objectName)
      throws StorageException {
    // define target blob object resource
    BlobInfo blobInfo = BlobInfo.newBuilder(BlobId.of(bucketName, objectName)).build();
    String fileName =
        objectName.contains("/")
            ? objectName.substring(objectName.lastIndexOf("/") + 1)
            : objectName;

    // Add response-content-disposition to force download with a specified filename)
    Map<String, String> extensionHeaders = new HashMap<>();
    extensionHeaders.put(
        "response-content-disposition", "attachment; filename=\"" + fileName + "\"");

    // generate signed URL
    URL url =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                gcsClient
                    .getStorageServiceWithProject(projectId)
                    .signUrl(
                        blobInfo,
                        gcsConfiguration.signedUrlGetDurationHours(),
                        TimeUnit.HOURS,
                        Storage.SignUrlOption.httpMethod(HttpMethod.GET),
                        Storage.SignUrlOption.withQueryParams(extensionHeaders),
                        Storage.SignUrlOption.withV4Signature()));

    String cleanSignedUrlString = cleanSignedUrl(url);
    logger.info("Generated GET signed URL: {}", cleanSignedUrlString);

    return url;
  }

  /**
   * Redact the X-Google-Signature element's value from the signed url and return the cleaned result
   * as a string.
   *
   * @param signedUrl
   * @return
   */
  public static String cleanSignedUrl(URL signedUrl) {
    String signedUrlString = signedUrl.toString();
    String[] signedUrlParts = signedUrlString.split("\\?");
    String elementDelimiter = "&";
    List<String> signedUrlElements = List.of(signedUrlParts[1].split(elementDelimiter));

    String signatureKey = "X-Goog-Signature";

    String cleanUrl = signedUrlString;
    if (signedUrlString.contains(signatureKey)) {
      String urlElementsWithoutSignature =
          signedUrlElements.stream()
              .filter(signedUrlElement -> !signedUrlElement.contains(signatureKey))
              .collect(Collectors.joining(elementDelimiter));
      cleanUrl =
          signedUrlParts[0]
              + "?"
              + urlElementsWithoutSignature
              + elementDelimiter
              + signatureKey
              + "=REDACTED";
    }
    return cleanUrl;
  }

  interface GcsAction<T> {
    T execute();
  }

  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, GcsAction<T> action) {
    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (StorageException e) {
            // Note: GCS' StorageException contains retryable exceptions - not sure how to handle
            throw new GcsServiceException("Error executing GCS action", e);
          }
        });
  }
}
