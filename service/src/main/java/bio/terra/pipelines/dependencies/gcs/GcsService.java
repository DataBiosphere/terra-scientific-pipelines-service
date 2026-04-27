package bio.terra.pipelines.dependencies.gcs;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.external.GcsConfiguration;
import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.service.exception.RequesterPaysBucketException;
import com.google.cloud.storage.*;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.charset.StandardCharsets;
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
  private static final String BUCKET_RESOURCE_FORMAT_FOR_LOGS = "bucket %s";

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
   * Helper method to retrieve a Blob object for a given GCS file.
   *
   * @param gcsFile GcsFile object representing the GCS path to the file
   * @param accessToken for the calling user. Pass null to use application default credentials.
   * @return Blob object if the file exists and is accessible, null otherwise
   * @throws RequesterPaysBucketException if the bucket is a requester pays bucket
   */
  public Blob getFileBlob(GcsFile gcsFile, String accessToken) {
    BlobId blobId = BlobId.fromGsUtilUri(gcsFile.getFullPath());
    try {
      // attempt to retrieve a client-side representation of the blob with minimal metadata
      return gcsClient
          .getStorageService(accessToken)
          .get(
              blobId, Storage.BlobGetOption.fields(Storage.BlobField.NAME, Storage.BlobField.SIZE));
    } catch (StorageException e) {
      if (e.getMessage().contains("Bucket is a requester pays bucket")) {
        throw new RequesterPaysBucketException(
            "The bucket for file '%s' is a requester pays bucket, which is not currently supported. Please turn off requester pays or submit data from a non-requester pays bucket."
                .formatted(gcsFile.getFullPath()));
      } else {
        logger.error(
            "An error occurred retrieving GCS file metadata for path `{}`. Error: {}",
            gcsFile.getFullPath(),
            e.getMessage());
        return null;
      }
    }
  }

  /** Get a GCS Blob object for a given GcsFile, with retries. Returns null if no blob found. */
  private Blob getBlob(GcsFile gcsFile) {
    BlobId blobId = BlobId.of(gcsFile.getBucketName(), gcsFile.getObjectName());

    return executionWithRetryTemplate(
        listenerResetRetryTemplate, () -> gcsClient.getStorageService().get(blobId));
  }

  /**
   * Get a BufferedReader for reading the contents of a GCS file. Note: the caller is responsible
   * for closing the BufferedReader after use to free resources. This method assumes the file
   * content is text-based and encoded in UTF-8.
   *
   * @param gcsFile GcsFile object representing the GCS file to read
   * @return
   */
  public BufferedReader getBufferedReaderForGcsTextFile(GcsFile gcsFile) {
    Blob blob = getBlob(gcsFile);
    if (blob == null) {
      logger.error("Blob not found for GcsFile: {}", gcsFile.getFullPath());
      throw new GcsServiceException("GCS file not found");
    }
    BufferedInputStream bis = new BufferedInputStream(Channels.newInputStream(blob.reader()));
    // bridge from byte streams (BufferedInputStream) to character streams (InputStreamReader)
    InputStreamReader isr = new InputStreamReader(bis, StandardCharsets.UTF_8);
    // wrap the InputStreamReader in a BufferedReader to use readLine()
    return new BufferedReader(isr);
  }

  /**
   * Helper method to retrieve the size of a file in GCS in bytes.
   *
   * @param gcsFilePath the full GCS path to the file (e.g. gs://my-bucket/path/to/file.txt)
   * @return the size of the file in bytes
   * @throws InternalServerErrorException if the file does not exist
   * @throws RequesterPaysBucketException if the bucket is requester pays
   */
  public Long getFileSizeInBytes(String gcsFilePath) {
    GcsFile gcsFile = new GcsFile(gcsFilePath);
    Blob fileBlob = getFileBlob(gcsFile, null);

    if (fileBlob == null) {
      throw new InternalServerErrorException(
          "An error occurred while retrieving file size for '%s'.".formatted(gcsFilePath));
    }

    Long size = fileBlob.getSize();
    logger.debug("Retrieved file size for '{}': {} bytes", gcsFile.getFileName(), size);
    return size;
  }

  /**
   * Check if a given user has read access to a GCS file.
   *
   * @param gcsFile GcsFile object representing the GCS path to the file
   * @param accessToken for the calling user. pass null to use application default credentials.
   * @return boolean whether the caller has read access to the GCS file
   * @throws RequesterPaysBucketException if the bucket is a requester pays bucket
   */
  private boolean hasFileReadAccess(GcsFile gcsFile, String accessToken) {
    try {
      Blob blob = getFileBlob(gcsFile, accessToken);
      if (blob != null && blob.exists()) {
        return true;
      }
    } catch (StorageException e) {
      logger.error("An error occurred checking blob access: {}", e.getMessage());
      return false;
    }
    return false;
  }

  public boolean serviceHasBucketReadAccess(String bucketName) {
    boolean hasAccess = hasBucketReadAccess(bucketName, null);
    logAccessCheckResult(
        SERVICE_SUBJECT_FOR_LOGS,
        "read",
        BUCKET_RESOURCE_FORMAT_FOR_LOGS.formatted(bucketName),
        hasAccess);
    return hasAccess;
  }

  public boolean userHasBucketReadAccess(String bucketName, @NonNull String accessToken) {
    boolean hasAccess = hasBucketReadAccess(bucketName, accessToken);
    logAccessCheckResult(
        USER_SUBJECT_FOR_LOGS,
        "read",
        BUCKET_RESOURCE_FORMAT_FOR_LOGS.formatted(bucketName),
        hasAccess);
    return hasAccess;
  }

  private boolean hasBucketReadAccess(String bucketName, String accessToken) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> {
          try {
            List<Boolean> accessResult =
                gcsClient
                    .getStorageService(accessToken)
                    .testIamPermissions(bucketName, List.of("storage.objects.get"));
            return accessResult.size() == 1 && accessResult.get(0);
          } catch (StorageException e) {
            if (e.getMessage().contains("Bucket is a requester pays bucket")) {
              throw new RequesterPaysBucketException(
                  "The bucket '%s' is a requester pays bucket, which is not currently supported. Please turn off requester pays or submit data from a non-requester pays bucket."
                      .formatted(bucketName));
            } else {
              throw e;
            }
          }
        });
  }

  public boolean serviceHasBucketWriteAccess(String bucketName) {
    boolean hasAccess = hasBucketWriteAccess(bucketName, null);
    logAccessCheckResult(
        SERVICE_SUBJECT_FOR_LOGS,
        "write",
        BUCKET_RESOURCE_FORMAT_FOR_LOGS.formatted(bucketName),
        hasAccess);
    return hasAccess;
  }

  public boolean userHasBucketWriteAccess(String bucketName, @NonNull String accessToken) {
    boolean hasAccess = hasBucketWriteAccess(bucketName, accessToken);
    logAccessCheckResult(
        USER_SUBJECT_FOR_LOGS,
        "write",
        BUCKET_RESOURCE_FORMAT_FOR_LOGS.formatted(bucketName),
        hasAccess);
    return hasAccess;
  }

  private boolean hasBucketWriteAccess(String bucketName, String accessToken) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> {
          List<Boolean> accessResult =
              gcsClient
                  .getStorageService(accessToken)
                  .testIamPermissions(
                      bucketName, List.of("storage.objects.create", "storage.objects.delete"));
          return accessResult.size() == 2 && accessResult.stream().allMatch(Boolean::booleanValue);
        });
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
    BlobId sourceBlob = BlobId.fromGsUtilUri(sourceUri.getFullPath());
    BlobId targetBlob = BlobId.fromGsUtilUri(targetUri.getFullPath());

    executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> {
          Storage storage = gcsClient.getStorageService();
          storage.copy(
              Storage.CopyRequest.newBuilder().setSource(sourceBlob).setTarget(targetBlob).build());

          logger.info(
              "Copied object from {} to {}", sourceUri.getFullPath(), targetUri.getFullPath());
          return targetBlob;
        });
  }

  /**
   * Delete a GCS object at the specified location.
   *
   * @param targetUri the full gs:// uri of the object to delete
   */
  public void deleteObject(GcsFile targetUri) throws StorageException {
    BlobId targetBlob = BlobId.fromGsUtilUri(targetUri.getFullPath());

    // The GCS delete method returns false if the object was not found, true if it was deleted.
    // If the call fails for some other reason, it will throw an exception which will be retried.
    boolean deleted =
        executionWithRetryTemplate(
            listenerResetRetryTemplate, () -> gcsClient.getStorageService().delete(targetBlob));
    if (deleted) {
      logger.info("Deleted object {}", targetUri.getFullPath());
    } else {
      logger.warn(
          "Object {} was not found for deletion. It may have already been deleted.",
          targetUri.getFullPath());
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
