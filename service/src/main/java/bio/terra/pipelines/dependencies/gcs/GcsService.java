package bio.terra.pipelines.dependencies.gcs;

import bio.terra.pipelines.app.configuration.external.GcsConfiguration;
import com.google.cloud.storage.BlobId;
import com.google.cloud.storage.BlobInfo;
import com.google.cloud.storage.HttpMethod;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/** class to encapsulate interacting with GCS client */
@Service
public class GcsService {

  private final GcsClient gcsClient;
  private final GcsConfiguration gcsConfiguration;

  private static final Logger logger = LoggerFactory.getLogger(GcsService.class);

  public GcsService(GcsClient gcsClient, GcsConfiguration gcsConfiguration) {
    this.gcsClient = gcsClient;
    this.gcsConfiguration = gcsConfiguration;
  }

  /**
   * Generates and returns a PUT (write-only) signed url for a specific object in a bucket.
   * See documentation on signed urls <a href="https://cloud.google.com/storage/docs/access-control/signed-urls">here</a>.
   *
   * The output URL can be used with a curl command to upload an object to the destination:
   * `curl -X PUT -H 'Content-Type: application/octet-stream' --upload-file my-file '{url}'`
   *
   * @param projectId Google project id
   * @param bucketName without a prefix
   * @param objectName should include the full path of the object (subdirectories + file name)
   * @return url that can be used to write an object to GCS
   */
  public URL generatePutObjectSignedUrl(String projectId, String bucketName, String objectName)
      throws StorageException {
    Storage storageService = gcsClient.getStorageService(projectId);

    // define target blob object resource
    BlobInfo blobInfo = BlobInfo.newBuilder(BlobId.of(bucketName, objectName)).build();

    // generate signed URL
    Map<String, String> extensionHeaders = new HashMap<>();
    extensionHeaders.put("Content-Type", "application/octet-stream");

    URL url =
        storageService.signUrl(
            blobInfo,
            gcsConfiguration.signedUrlPutDurationHours(),
            TimeUnit.HOURS,
            Storage.SignUrlOption.httpMethod(HttpMethod.PUT),
            Storage.SignUrlOption.withExtHeaders(extensionHeaders),
            Storage.SignUrlOption.withV4Signature());

    logger.info("Generated PUT signed URL:%s".formatted(url.toString()));

    return url;
  }
}
