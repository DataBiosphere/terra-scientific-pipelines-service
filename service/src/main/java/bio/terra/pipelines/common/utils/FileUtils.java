package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.InternalServerErrorException;
import java.util.UUID;

/** A collection of utilities and constants useful for files. */
public class FileUtils {
  private FileUtils() {
    throw new IllegalStateException("Utility class");
  }

  private static final String USER_PROVIDED_FILE_INPUT_DIRECTORY = "user-input-files";

  /**
   * Extract the blob name from the full file path, using the workspaceId as a delimiter.
   *
   * <p>For example, `https://lz123.blob.core.windows.net/sc-{workspaceId}/path/to/file` becomes
   * `path/to/file`
   *
   * @param blobHttpUrl
   * @param workspaceId
   * @return blobName
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
   * Construct the destination blob name for a user-provided file input.
   *
   * <p>For example, file `local/path/to/file.txt` for jobId `1234` returns
   * `user-input-files/1234/file.txt`
   *
   * @param jobId
   * @param userProvidedFileInputValue
   * @return blobName
   */
  public static String constructDestinationBlobNameForUserInputFile(
      UUID jobId, String userProvidedFileInputValue) {
    String userProvidedFileName = getFileNameFromFullPath(userProvidedFileInputValue);
    return "%s/%s/%s".formatted(USER_PROVIDED_FILE_INPUT_DIRECTORY, jobId, userProvidedFileName);
  }

  /**
   * Extract the storage container url from a full path or SAS url, using the workspaceId as a
   * delimiter.
   *
   * <p>For example, `https://lz123.blob.core.windows.net/sc-{workspaceId}/path/to/file` becomes
   * `https://lz123.blob.core.windows.net/sc-{workspaceId}`
   *
   * @param sasUrl
   * @param workspaceId
   * @return baseStorageUrl
   */
  public static String getStorageContainerUrlFromSasUrl(String sasUrl, UUID workspaceId) {
    if (!sasUrl.contains(workspaceId.toString())) {
      throw new InternalServerErrorException(
          "File path and workspaceId do not match. Cannot extract base storage url.");
    }
    return sasUrl.substring(
        0, sasUrl.indexOf(workspaceId.toString()) + workspaceId.toString().length());
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
