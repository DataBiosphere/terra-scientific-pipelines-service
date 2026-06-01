package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class PipelinesEnumTest extends BaseTest {

  @Test
  void buildPipelineKey() {
    assertEquals(
        "array_imputation_v1", PipelinesEnum.buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, 1));
    assertEquals(
        "low_pass_imputation_v3",
        PipelinesEnum.buildPipelineKey(PipelinesEnum.LOW_PASS_IMPUTATION, 3));
  }

  @Test
  void nameFromPipelineKeyValid() {
    assertEquals(
        PipelinesEnum.ARRAY_IMPUTATION, PipelinesEnum.nameFromPipelineKey("array_imputation_v2"));
    assertEquals(
        PipelinesEnum.LOW_PASS_IMPUTATION,
        PipelinesEnum.nameFromPipelineKey("low_pass_imputation_v1"));
  }

  @Test
  void nameFromPipelineKeyThrowsNotFound() {
    assertThrows(
        NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey("ARRAY_IMPUTATION_v1"));
    assertThrows(
        NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey("ARRAY_IMPUTATION_V1"));
    assertThrows(
        NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey("array_imputation_V1"));
    assertThrows(
        NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey("array_imputation_v"));
    assertThrows(
        NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey("unknown_pipeline_v1"));
    assertThrows(NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey("nounderscorev"));
    assertThrows(NotFoundException.class, () -> PipelinesEnum.nameFromPipelineKey(null));
  }

  @Test
  void versionFromPipelineKeyValid() {
    assertEquals(2, PipelinesEnum.versionFromPipelineKey("array_imputation_v2"));
    assertEquals(10, PipelinesEnum.versionFromPipelineKey("low_pass_imputation_v10"));
  }

  @Test
  void versionFromPipelineKeyThrowsNotFound() {
    assertThrows(
        NotFoundException.class,
        () -> PipelinesEnum.versionFromPipelineKey("array_imputation_vabc"));
    assertThrows(
        NotFoundException.class, () -> PipelinesEnum.versionFromPipelineKey("noprefixatall"));
    assertThrows(NotFoundException.class, () -> PipelinesEnum.versionFromPipelineKey(null));
  }

  @Test
  void buildAndParseRoundTrip() {
    String key = PipelinesEnum.buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, 5);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, PipelinesEnum.nameFromPipelineKey(key));
    assertEquals(5, PipelinesEnum.versionFromPipelineKey(key));
  }
}
