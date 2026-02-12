package bio.terra.pipelines.common;

import com.google.cloud.storage.BlobId;
import java.util.Objects;
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

  // we override equals() and hashCode() so that we can compare GcsFile objects in tests; fullPath
  // is the only unique identifier for a GcsFile, so we use that for equality and hash code
  // generation
  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    GcsFile gcsFile = (GcsFile) o;
    return Objects.equals(fullPath, gcsFile.fullPath);
  }

  @Override
  public int hashCode() {
    return Objects.hash(fullPath);
  }
}
