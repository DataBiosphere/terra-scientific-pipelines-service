package bio.terra.pipelines.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

public class GcsFileTest extends BaseTest {
  @Test
  public void gcsFile() {
    String uri = "gs://my-bucket/path/to/file.txt";
    GcsFile gcsFile = new GcsFile(uri);
    assertEquals("my-bucket", gcsFile.getBucketName());
    assertEquals("path/to/file.txt", gcsFile.getObjectName());
    assertEquals("file.txt", gcsFile.getFileName());
    assertEquals(uri, gcsFile.getFullPath());

    assertThrows(IllegalArgumentException.class, () -> new GcsFile("not-a-uri"));
    assertThrows(IllegalArgumentException.class, () -> new GcsFile("gs://this-is-just-a-bucket"));
    assertThrows(
        IllegalArgumentException.class, () -> new GcsFile("gs://this-is-also-just-a-bucket/"));
    assertThrows(IllegalArgumentException.class, () -> new GcsFile("s3://aws-bucket/file"));
  }
}
