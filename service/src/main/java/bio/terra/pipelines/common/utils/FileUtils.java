package bio.terra.pipelines.common.utils;

import java.net.URI;
import java.net.URISyntaxException;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** A collection of utilities and constants useful for files. */
public class FileUtils {
  private FileUtils() {
    throw new IllegalStateException("Utility class");
  }

  private static final String USER_PROVIDED_FILE_INPUT_DIRECTORY = "user-input-files";
  public static final String GCP_STORAGE_PROTOCOL = "gs://";
  private static final Pattern GCS_BUCKET_PATTERN =
      //      Pattern.compile("^gs://([a-z0-9][a-z0-9._-]{1,61}[a-z0-9])/(.+)$");
      Pattern.compile("^gs://([^/]+)/.+$");

  /**
   * Extracts the bucket name from a GCS path. For example, `gs://my-bucket/path/to/file.txt`
   * returns `my-bucket`. Returns null if the bucket name cannot be extracted.
   */
  public static String extractGcsBucketName(String item) {
    Matcher m = GCS_BUCKET_PATTERN.matcher(item);
    return m.matches() ? m.group(1) : null;
  }

  /**
   * Construct the destination blob name for a local user-provided file input.
   *
   * <p>For example, file `local/path/to/file.txt` for jobId `1234` returns
   * `user-input-files/1234/file.txt`
   *
   * @param jobId the jobID for the pipelineRun
   * @param userProvidedFileInputValue the file input value provided by the user, which is expected
   *     to be a local file path
   * @return blobName
   */
  public static String constructDestinationBlobNameForUserInputFile(
      UUID jobId, String userProvidedFileInputValue) {
    String userProvidedFileName = getFileNameFromFullPath(userProvidedFileInputValue);
    return constructFilePath(
        constructFilePath(USER_PROVIDED_FILE_INPUT_DIRECTORY, jobId.toString()),
        userProvidedFileName);
  }

  /**
   * Construct the full GCS file path for a local user-provided file input, which includes the
   * bucket and the custom blob name/path.
   *
   * @param bucketName can include the gs:// protocol or not
   * @param jobId the jobID for the pipelineRun
   * @param userProvidedFileInputValue the file input value provided by the user, which is expected
   *     to be a local file path
   * @return String fully qualified gsutil uri pointing to the file in the workspace
   */
  public static String constructGcsFilePathForUserLocalInputFile(
      String bucketName, UUID jobId, String userProvidedFileInputValue) {
    if (!bucketName.startsWith(GCP_STORAGE_PROTOCOL)) {
      bucketName = GCP_STORAGE_PROTOCOL + bucketName;
    }
    String blobName =
        constructDestinationBlobNameForUserInputFile(jobId, userProvidedFileInputValue);
    return constructFilePath(bucketName, blobName);
  }

  /** Determine the file location type from a file path. */
  public static FileLocationTypeEnum getFileLocationType(String filePath) {
    try {
      URI uri = new URI(filePath);
      String scheme = uri.getScheme();
      if (scheme == null) {
        return FileLocationTypeEnum.LOCAL;
      } else if (scheme.equals("gs")) {
        return FileLocationTypeEnum.GCS;
      } else {
        return FileLocationTypeEnum.UNSUPPORTED;
      }
    } catch (
        URISyntaxException
            e) { // catches strings with characters like spaces, these should not have passed
      // validation anyway
      return FileLocationTypeEnum.UNSUPPORTED;
    }
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
