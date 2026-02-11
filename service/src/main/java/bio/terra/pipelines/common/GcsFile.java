package bio.terra.pipelines.common;

import com.google.cloud.storage.BlobId;
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class GcsFile {
  private final String bucketName;
  private final String objectName;
  private final String fullPath;
  private final String fileName;

  public GcsFile(String fullPath) {
    this.fullPath = fullPath;
    BlobId blobId = BlobId.fromGsUtilUri(fullPath);
    this.bucketName = blobId.getBucket();
    this.objectName = blobId.getName();
    this.fileName =
        objectName.contains("/")
            ? objectName.substring(objectName.lastIndexOf("/") + 1)
            : objectName;
  }
}
