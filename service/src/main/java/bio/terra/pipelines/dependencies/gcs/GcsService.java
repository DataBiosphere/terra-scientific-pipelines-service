package bio.terra.pipelines.dependencies.gcs;

import bio.terra.pipelines.app.configuration.external.GcsConfiguration;
import bio.terra.pipelines.common.GcsFile;
import com.google.cloud.storage.*;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import lombok.NonNull;
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

  private static final String SERVICE_SUBJECT_FOR_LOGS = "Teaspoons service account";
  private static final String USER_SUBJECT_FOR_LOGS = "User";
  private static final String FILE_RESOURCE_FORMAT_FOR_LOGS = "file %s";

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

  public boolean serviceHasFileReadAccess(GcsFile gcsFile) {
    boolean hasAccess = hasFileReadAccess(gcsFile, null);
    logAccessCheckResult(
        SERVICE_SUBJECT_FOR_LOGS,
        "read",
        FILE_RESOURCE_FORMAT_FOR_LOGS.formatted(gcsFile.getFullPath()),
        hasAccess);
    return hasAccess;
  }

  public boolean userHasFileReadAccess(GcsFile gcsFile, @NonNull String accessToken) {
    boolean hasAccess = hasFileReadAccess(gcsFile, accessToken);
    logAccessCheckResult(
        USER_SUBJECT_FOR_LOGS,
        "read",
        FILE_RESOURCE_FORMAT_FOR_LOGS.formatted(gcsFile.getFullPath()),
        hasAccess);
    return hasAccess;
  }

  /**
   * Check if a given user has read access to a GCS file.
   *
   * @param gcsFile GcsFile object representing the GCS path to the file
   * @param accessToken for the calling user. pass null to use application default credentials.
   * @return boolean whether the caller has read access to the GCS file
   */
  private boolean hasFileReadAccess(GcsFile gcsFile, String accessToken) {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
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
      logger.error("An error occurred checking blob access: {}", e.getMessage());
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
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object (subdirectories + file name)
   * @return url that can be used to write an object to GCS
   */
  public URL generatePutObjectSignedUrl(String bucketName, String objectName)
      throws StorageException {
    return generateUploadSignedUrl(bucketName, objectName, false);
  }

  /**
   * Generates and returns a resumable POST (write-only) signed url for a specific object in a
   * bucket. See documentation on signed urls <a
   * href="https://cloud.google.com/storage/docs/access-control/signed-urls">here</a> and resumable
   * uploads <a
   * href="https://cloud.google.com/storage/docs/performing-resumable-uploads#initiate-session">here</a>.
   *
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object (subdirectories + file name)
   * @return url that can be used to initiate a resumable object upload to GCS
   */
  public URL generateResumablePostObjectSignedUrl(String bucketName, String objectName)
      throws StorageException {
    return generateUploadSignedUrl(bucketName, objectName, true);
  }

  private URL generateUploadSignedUrl(String bucketName, String objectName, boolean isResumable)
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
                    .getStorageService()
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
   * @param gcsFile object representing the GCS file for which to generate a signed URL
   * @return url that can be used to download an object to GCS
   */
  public URL generateGetObjectSignedUrl(GcsFile gcsFile) throws StorageException {
    // define target blob object resource
    BlobInfo blobInfo = BlobInfo.newBuilder(BlobId.fromGsUtilUri(gcsFile.getFullPath())).build();

    // Add response-content-disposition to force download with a specified filename)
    Map<String, String> extensionHeaders = new HashMap<>();
    extensionHeaders.put(
        "response-content-disposition", "attachment; filename=\"" + gcsFile.getFileName() + "\"");

    // generate signed URL
    URL url =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                gcsClient
                    .getStorageService()
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
   * Copy a GCS object from one location to another.
   *
   * @param sourceUri full GCS path of the source file
   * @param targetUri full GCS path of the target file
   */
  public void copyObject(GcsFile sourceUri, GcsFile targetUri) throws StorageException {

    BlobId source = BlobId.fromGsUtilUri(sourceUri.getFullPath());
    BlobId target = BlobId.fromGsUtilUri(targetUri.getFullPath());

    executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> {
          Storage storage = gcsClient.getStorageService();
          storage.copy(
              Storage.CopyRequest.newBuilder().setSource(source).setTarget(target).build());
          Blob copiedObject = storage.get(target);

          logger.info(
              "Copied object {} from bucket {} to {} in bucket {}",
              source.getName(),
              source.getBucket(),
              target.getName(),
              copiedObject.getBucket());
          return target;
        });
  }

  /**
   * Delete a GCS object at the specified location.
   *
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object
   */
  public void deleteObject(String bucketName, String objectName) throws StorageException {
    BlobId blobId = BlobId.of(bucketName, objectName);
    boolean deleted =
        executionWithRetryTemplate(
            listenerResetRetryTemplate, () -> gcsClient.getStorageService().delete(blobId));
    if (deleted) {
      logger.info("Deleted object {} in bucket {}", objectName, bucketName);
    } else {
      logger.warn(
          "Object {} in bucket {} was not found for deletion. It may have already been deleted.",
          objectName,
          bucketName);
    }
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
