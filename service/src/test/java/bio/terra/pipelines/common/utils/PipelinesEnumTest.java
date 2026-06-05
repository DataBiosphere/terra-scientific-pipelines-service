package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class PipelinesEnumTest extends BaseTest {
  @Test
  void enumFromConfigKey() {
    assertEquals(
        PipelinesEnum.ARRAY_IMPUTATION, PipelinesEnum.enumFromConfigKeyValue("arrayImputation"));
    assertEquals(
        PipelinesEnum.LOW_PASS_IMPUTATION,
        PipelinesEnum.enumFromConfigKeyValue("lowPassImputation"));
  }
}
