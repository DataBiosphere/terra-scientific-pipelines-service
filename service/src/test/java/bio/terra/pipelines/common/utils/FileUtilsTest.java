package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.common.exception.BadRequestException;
import bio.terra.pipelines.testutils.BaseTest;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.UUID;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

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

  private static Stream<Arguments> parseTsvInputs() {
    return Stream.of(
        // arguments: string to be converted to an InputStream, expected output (null if expect
        // error)
        arguments(
            "one\ttwo\tthree\nfour\tfive\tsix",
            List.of(new String[] {"one", "two", "three"}, new String[] {"four", "five", "six"})),
        arguments(
            "one\ttwo\tthree\nfour\tfive\tsix\n",
            List.of(new String[] {"one", "two", "three"}, new String[] {"four", "five", "six"})),
        arguments("", List.of()),
        arguments( // allow missing values
            "one\ttwo\tthree\n\t\tfour",
            List.of(new String[] {"one", "two", "three"}, new String[] {"", "", "four"})),
        arguments(
            "just a lot of text", Collections.singletonList(new String[] {"just a lot of text"})),
        arguments(
            "col1,col2,val1,val2", Collections.singletonList(new String[] {"col1,col2,val1,val2"})),
        // failure cases
        arguments("col1\tcol2\tcol3\nval1\tval2", null));
  }

  @ParameterizedTest
  @MethodSource("parseTsvInputs")
  void parseTsv(String inputStreamString, List<String[]> expectedOutput) throws IOException {
    try (InputStream inputStream = createInputStreamForTesting(inputStreamString)) {
      if (expectedOutput != null) {
        List<String[]> actualOutput = FileUtils.parseTsv(inputStream);
        System.out.println("Actual Output: " + Arrays.deepToString(actualOutput.toArray()));
        System.out.println("Expected Output: " + Arrays.deepToString(expectedOutput.toArray()));

        assertEquals(expectedOutput.size(), actualOutput.size());
        for (int i = 0; i < expectedOutput.size(); i++) {
          assertArrayEquals(expectedOutput.get(i), actualOutput.get(i));
        }
      } else {
        assertThrows(BadRequestException.class, () -> FileUtils.parseTsv(inputStream));
      }
    }
  }

  private static Stream<Arguments> getItemsFromManifestLinesInputs() {
    return Stream.of(
        // arguments: input List of String arrays, expected output List of Strings
        arguments(
            List.of(new String[] {"col1", "col2", "col3"}, new String[] {"val1", "val2", "val3"}),
            List.of("col1", "col2", "col3", "val1", "val2", "val3")),
        arguments(List.of(), List.of()),
        arguments(
            List.of(new String[] {"val1"}, new String[] {"val2", "val3"}),
            List.of("val1", "val2", "val3")));
  }

  @ParameterizedTest
  @MethodSource("getItemsFromManifestLinesInputs")
  void getItemsFromManifestLines(List<String[]> manifestLines, List<String> expectedOutput) {
    assertEquals(expectedOutput, FileUtils.getItemsFromManifestLines(manifestLines));
  }

  /** Helper method to create an InputStream from a string for testing purposes. */
  public static InputStream createInputStreamForTesting(String testData) {
    byte[] bytes = testData.getBytes(StandardCharsets.UTF_8);
    return new ByteArrayInputStream(bytes);
  }
}
