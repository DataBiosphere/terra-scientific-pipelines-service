package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.testutils.BaseTest;
import java.util.UUID;
import org.junit.jupiter.api.Test;

class FileUtilsTest extends BaseTest {

  @Test
  void getBlobNameFromTerraWorkspaceStorageHttpUrl() {
    String fullPath =
        "https://lze96253b07f13c61ef712bb.blob.core.windows.net/sc-68a43bd8-e744-4f1e-87a5-c44ecef157a3/workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    UUID controlWorkspaceId = UUID.fromString("68a43bd8-e744-4f1e-87a5-c44ecef157a3");
    String expectedBlobName =
        "workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    assertEquals(
        expectedBlobName,
        FileUtils.getBlobNameFromTerraWorkspaceStorageHttpUrl(fullPath, controlWorkspaceId));
  }

  @Test
  void getBlobNameFromTerraWorkspaceStorageHttpUrlDifferentWorkspace() {
    String fullPath =
        "https://lze96253b07f13c61ef712bb.blob.core.windows.net/sc-68a43bd8-e744-4f1e-87a5-c44ecef157a3/workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    UUID wrongWorkspaceId = UUID.fromString("11111111-1111-1111-1111-111111111111");

    assertThrows(
        InternalServerErrorException.class,
        () -> FileUtils.getBlobNameFromTerraWorkspaceStorageHttpUrl(fullPath, wrongWorkspaceId));
  }

  @Test
  void constructBlobNameForUserInputFile() {
    UUID jobId = UUID.randomUUID();
    String userProvidedFileInputValue = "local/path/to/file.txt";
    String expectedBlobName = "user-input-files/%s/file.txt".formatted(jobId);
    assertEquals(
        expectedBlobName,
        FileUtils.constructBlobNameForUserInputFile(jobId, userProvidedFileInputValue));
  }

  @Test
  void getBaseStorageUrlFromSasUrl() {
    String sasUrl =
        "https://lz123.blob.core.windows.net/sc-68a43bd8-e744-4f1e-87a5-c44ecef157a3/workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    UUID controlWorkspaceId = UUID.fromString("68a43bd8-e744-4f1e-87a5-c44ecef157a3");
    String expectedBaseStorageUrl =
        "https://lz123.blob.core.windows.net/sc-68a43bd8-e744-4f1e-87a5-c44ecef157a3";
    assertEquals(
        expectedBaseStorageUrl, FileUtils.getBaseStorageUrlFromSasUrl(sasUrl, controlWorkspaceId));
  }

  @Test
  void getBaseStorageUrlFromSasUrlDifferentWorkspace() {
    String sasUrl =
        "https://lz123.blob.core.windows.net/sc-68a43bd8-e744-4f1e-87a5-c44ecef157a3/workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    UUID wrongWorkspaceId = UUID.fromString("11111111-1111-1111-1111-111111111111");

    assertThrows(
        InternalServerErrorException.class,
        () -> FileUtils.getBaseStorageUrlFromSasUrl(sasUrl, wrongWorkspaceId));
  }

  @Test
  void getFileNameFromFullPath() {
    String fullPath = "path/to/file.txt";
    String expectedFileName = "file.txt";
    assertEquals(expectedFileName, FileUtils.getFileNameFromFullPath(fullPath));
  }
}
