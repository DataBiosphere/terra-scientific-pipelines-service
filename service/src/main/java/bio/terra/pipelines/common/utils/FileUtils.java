package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.InternalServerErrorException;
import java.util.UUID;

/** A collection of utilities and constants useful for files. */
public class FileUtils {
  private FileUtils() {
    throw new IllegalStateException("Utility class");
  }

  /**
   * Extract the blob name from the full file path, using the workspaceId as a delimiter.
   *
   * <p>For example, `https://lz123.blob.core.windows.net/sc-{workspaceId}/path/to/file` becomes
   * `path/to/file`
   */
  public static String getBlobNameFromTerraWorkspaceStorageHttpUrl(
      String blobHttpUrl, UUID workspaceId) {
    if (!blobHttpUrl.contains(workspaceId.toString())) {
      throw new InternalServerErrorException(
          "File path and workspaceId do not match. Cannot extract blob name.");
    }
    return blobHttpUrl.substring(
        blobHttpUrl.indexOf(workspaceId.toString()) + workspaceId.toString().length() + 1);
  }

  /**
   * Extract the file name from a full path.
   *
   * <p>For example, `path/to/file.txt` becomes `file.txt`
   *
   * @param fullPath the full path to the file
   * @return fileName
   */
  public static String getFileNameFromFullPath(String fullPath) {
    return fullPath.substring(fullPath.lastIndexOf('/') + 1);
  }
}
