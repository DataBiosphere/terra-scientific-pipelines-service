package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.enumFromPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.versionFromPipelineKey;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class PipelineKeyUtilsTest extends BaseTest {

  @Test
  void buildPipelineKeyTest() {
    assertEquals("array_imputation_v1", buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, 1));
    assertEquals("low_pass_imputation_v3", buildPipelineKey(PipelinesEnum.LOW_PASS_IMPUTATION, 3));
  }

  @Test
  void enumFromPipelineKeyValid() {
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, enumFromPipelineKey("array_imputation_v2"));
    assertEquals(PipelinesEnum.LOW_PASS_IMPUTATION, enumFromPipelineKey("low_pass_imputation_v1"));
  }

  @Test
  void enumFromPipelineKeyThrowsNotFound() {
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey("ARRAY_IMPUTATION_v1"));
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey("ARRAY_IMPUTATION_V1"));
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey("array_imputation_V1"));
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey("array_imputation_v"));
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey("unknown_pipeline_v1"));
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey("nounderscorev"));
    assertThrows(NotFoundException.class, () -> enumFromPipelineKey(null));
  }

  @Test
  void versionFromPipelineKeyValid() {
    assertEquals(2, versionFromPipelineKey("array_imputation_v2"));
    assertEquals(10, versionFromPipelineKey("low_pass_imputation_v10"));
  }

  @Test
  void versionFromPipelineKeyThrowsNotFound() {
    assertThrows(NotFoundException.class, () -> versionFromPipelineKey("array_imputation_vabc"));
    assertThrows(NotFoundException.class, () -> versionFromPipelineKey("noprefixatall"));
    assertThrows(NotFoundException.class, () -> versionFromPipelineKey(null));
  }

  @Test
  void buildAndParseRoundTrip() {
    String key = buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, 5);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, enumFromPipelineKey(key));
    assertEquals(5, versionFromPipelineKey(key));
  }
}
