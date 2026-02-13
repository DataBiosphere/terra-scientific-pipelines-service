package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.testutils.BaseTest;
import java.util.UUID;
import org.junit.jupiter.api.Test;

class FileUtilsTest extends BaseTest {

  @Test
  void constructBlobNameForUserInputFile() {
    UUID jobId = UUID.randomUUID();
    String userProvidedFileInputValue = "local/path/to/file.txt";
    String expectedBlobName = "user-input-files/%s/file.txt".formatted(jobId);
    assertEquals(
        expectedBlobName,
        FileUtils.constructDestinationBlobNameForUserInputFile(jobId, userProvidedFileInputValue));
  }

  @Test
  void getFileLocationType() {
    String gcsPath = "gs://bucket_name/path/to/file.txt";
    String awsPath = "s3://bucket_name/path/to/file.txt";
    String localPath = "/local/path/to/file.txt";
    String notAValidCloudPath = "htp://invalid/path/to/file.txt";
    String badCharacters = "gs://bucket_name/bad/pa th/to/file.txt";

    assertEquals(FileLocationTypeEnum.GCS, FileUtils.getFileLocationType(gcsPath));
    assertEquals(FileLocationTypeEnum.UNSUPPORTED, FileUtils.getFileLocationType(awsPath));
    assertEquals(FileLocationTypeEnum.LOCAL, FileUtils.getFileLocationType(localPath));
    assertEquals(
        FileLocationTypeEnum.UNSUPPORTED, FileUtils.getFileLocationType(notAValidCloudPath));
    assertEquals(FileLocationTypeEnum.UNSUPPORTED, FileUtils.getFileLocationType(badCharacters));
  }

  @Test
  void getFileNameFromFullPath() {
    String fullPath = "path/to/file.txt";
    String expectedFileName = "file.txt";
    assertEquals(expectedFileName, FileUtils.getFileNameFromFullPath(fullPath));
  }

  @Test
  void constructFilePath() {
    String pathWithoutEndingSlash = "path/to";
    String pathWithEndingSlash = "path/to/";
    String fileNameWithBeginningSlash = "/file.txt";
    String fileNameWithoutBeginningSlash = "file.txt";

    String expectedFullPath = "path/to/file.txt";

    // all combinations of paths and file names should produce expectedFullPath
    assertEquals(
        expectedFullPath,
        FileUtils.constructFilePath(pathWithoutEndingSlash, fileNameWithBeginningSlash));
    assertEquals(
        expectedFullPath,
        FileUtils.constructFilePath(pathWithEndingSlash, fileNameWithBeginningSlash));
    assertEquals(
        expectedFullPath,
        FileUtils.constructFilePath(pathWithoutEndingSlash, fileNameWithoutBeginningSlash));
    assertEquals(
        expectedFullPath,
        FileUtils.constructFilePath(pathWithEndingSlash, fileNameWithoutBeginningSlash));
  }
}
