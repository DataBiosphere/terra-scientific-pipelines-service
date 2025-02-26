package bio.terra.pipelines.db.entities;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class PipelineQuotaTest extends BaseTest {

  @Test
  void testPipelineQuotaConstructor() {
    PipelineQuota pipelineQuota =
        new PipelineQuota(PipelinesEnum.ARRAY_IMPUTATION, 10000, 500, QuotaUnitsEnum.SAMPLES);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipelineQuota.getPipelineName());
    assertEquals(10000, pipelineQuota.getDefaultQuota());
    assertEquals(500, pipelineQuota.getMinQuotaConsumed());
    assertEquals(QuotaUnitsEnum.SAMPLES, pipelineQuota.getQuotaUnits());
  }
}
