package bio.terra.pipelines.dependencies.gcs;

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

  private static final Logger logger = LoggerFactory.getLogger(GcsService.class);

  public GcsService(GcsClient gcsClient) {
    this.gcsClient = gcsClient;
  }

  /**
   * Signing a URL requires Credentials which implement ServiceAccountSigner. These can be set
   * explicitly using the Storage.SignUrlOption.signWith(ServiceAccountSigner) option. If you don't,
   * you could also pass a service account signer to StorageOptions, i.e.
   * StorageOptions().newBuilder().setCredentials(ServiceAccountSignerCredentials). In this example,
   * neither of these options are used, which means the following code only works when the
   * credentials are defined via the environment variable GOOGLE_APPLICATION_CREDENTIALS, and those
   * credentials are authorized to sign a URL. See the documentation for Storage.signUrl for more
   * details.
   */
  public URL generateV4PutObjectSignedUrl(String projectId, String bucketName, String objectName)
      throws StorageException {
    // String projectId = "my-project-id";
    // String bucketName = "my-bucket";
    // String objectName = "my-object";

    Storage storage = gcsClient.getStorageService(projectId);

    // Define Resource
    BlobInfo blobInfo = BlobInfo.newBuilder(BlobId.of(bucketName, objectName)).build();

    // Generate Signed URL
    Map<String, String> extensionHeaders = new HashMap<>();
    extensionHeaders.put("Content-Type", "application/octet-stream");

    URL url =
        storage.signUrl(
            blobInfo,
            15,
            TimeUnit.MINUTES,
            Storage.SignUrlOption.httpMethod(HttpMethod.PUT),
            Storage.SignUrlOption.withExtHeaders(extensionHeaders),
            Storage.SignUrlOption.withV4Signature());

    logger.info("Generated PUT signed URL:%s".formatted(url));

    return url;

    //        System.out.println("You can use this URL with any user agent, for example:");
    //        System.out.println(
    //                "curl -X PUT -H 'Content-Type: application/octet-stream' --upload-file my-file
    // '"
    //                        + url
    //                        + "'");
  }
}
