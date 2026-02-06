package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.InternalServerErrorException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.UUID;

/** A collection of utilities and constants useful for files. */
public class FileUtils {
  private FileUtils() {
    throw new IllegalStateException("Utility class");
  }

  private static final String USER_PROVIDED_FILE_INPUT_DIRECTORY = "user-input-files";

  /**
   * Extract the blob name from the full GCS file path, using a defined bucketName.
   *
   * <p>For example, with my-bucket as the bucketName, `gs://my-bucket/path/to/file` becomes
   * `path/to/file`.
   *
   * @param gcsUrl
   * @param bucketName
   * @return blobName
   */
  public static String getBlobNameFromGcsStorageUrl(String gcsUrl, String bucketName) {
    if (!gcsUrl.contains(bucketName)) {
      throw new InternalServerErrorException(
          "File path and bucketName do not match. Cannot extract blob name.");
    }
    return gcsUrl.substring(gcsUrl.indexOf(bucketName) + bucketName.length() + 1);
  }

  /**
   * Construct the destination blob name for a local user-provided file input.
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
   * Extract the GCS bucket from a full path
   *
   * <p>For example, `getBucketFromCloudPath("gs://bucket/path/to/file")` returns `"bucket"`
   *
   * @param cloudPath
   * @return bucketName or null if not a cloud path
   */
  public static String getBucketFromGcsCloudPath(String cloudPath) {
    if (getFileLocationType(cloudPath) == FileLocationTypeEnum.GCS) {
      return cloudPath.split("/")[2];
    }
    return null;
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
