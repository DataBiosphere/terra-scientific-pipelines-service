package bio.terra.pipelines.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class GcsFileTest extends BaseTest {
  @Test
  void gcsFile() {
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

  @Test
  void equals() {
    String path = "gs://bucket/path/to/file.txt";
    GcsFile first = new GcsFile(path);
    GcsFile sameAsFirst = new GcsFile(path);
    GcsFile different = new GcsFile("gs://other/path/to/file.txt");
    assertEquals(first, first);
    assertEquals(first, sameAsFirst);
    assertNotEquals(first, different);
    assertNotEquals(null, first);
    assertNotEquals(first, first.getFullPath());
  }

  @Test
  void hashCodeEquals() {
    String path = "gs://bucket/path/to/file.txt";
    GcsFile first = new GcsFile(path);
    GcsFile sameAsFirst = new GcsFile(path);
    GcsFile different = new GcsFile("gs://other/path/to/file.txt");
    assertEquals(first.hashCode(), first.hashCode());
    assertEquals(first.hashCode(), sameAsFirst.hashCode());
    assertNotEquals(first.hashCode(), different.hashCode());
  }
}
