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
   * Extract the blob name from the full GCP file path, using the workspaceStorageContainerName as a
   * delimiter.
   *
   * <p>For example, with workspaceStorageContainerName as the workspaceSubstringStart,
   * `gs://{workspaceDelimiter}/path/to/file` becomes `path/to/file`.
   *
   * @param blobUrl
   * @param workspaceStorageContainerName
   */
  public static String getBlobNameFromTerraWorkspaceStorageUrlGcp(
      String blobUrl, String workspaceStorageContainerName) {
    return getBlobNameFromTerraWorkspaceStorageUrl(blobUrl, workspaceStorageContainerName);
  }

  /**
   * Extract the blob name from the full file path, using a defined workspaceSubstringStart.
   *
   * @param blobUrl
   * @param workspaceSubstringStart
   * @return blobName
   */
  private static String getBlobNameFromTerraWorkspaceStorageUrl(
      String blobUrl, String workspaceSubstringStart) {
    if (!blobUrl.contains(workspaceSubstringStart)) {
      throw new InternalServerErrorException(
          "File path and workspaceSubstringStart do not match. Cannot extract blob name.");
    }
    return blobUrl.substring(
        blobUrl.indexOf(workspaceSubstringStart) + workspaceSubstringStart.length() + 1);
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
   * Determine whether a file path is a cloud path or local path. Return true if the path is in the
   * cloud (i.e. starts with "gs://"), false if it is a local path.
   */
  public static boolean isCloudFile(String filePath) {
    return filePath.startsWith("gs://");
  }

  /**
   * Extract the GCS bucket from a full path
   *
   * <p>For example, `getBucketFromCloudPath("gs://bucket/path/to/file")` returns `"bucket"`
   *
   * @param cloudPath
   * @return bucketName or null if not a cloud path
   */
  public static String getBucketFromCloudPath(String cloudPath) {
    if (isCloudFile(cloudPath)) {
      return cloudPath.split("/")[2];
    }
    return null;
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

  /**
   * Construct a path to a file given a base path and a file name. Adds a '/' between the base path
   * and the file name if the base path does not end with a '/' or the file name does not begin with
   * a '/'. Ensures there is only one '/' if both base and file paths contain a '/'.
   */
  public static String constructFilePath(String basePath, String fileName) {
    String pathDelimiter = "/";
    String basePathWithoutTrailingDelimiter =
        basePath.replaceAll("%s$".formatted(pathDelimiter), "");
    String fileNameWithoutLeadingDelimiter =
        fileName.replaceAll("^%s".formatted(pathDelimiter), "");
    return "%s/%s".formatted(basePathWithoutTrailingDelimiter, fileNameWithoutLeadingDelimiter);
  }
}
