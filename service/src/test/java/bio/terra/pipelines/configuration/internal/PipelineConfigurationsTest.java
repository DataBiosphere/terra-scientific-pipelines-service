package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.math.BigDecimal;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelineConfigurationsTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineConfigurations pipelineConfigurations;
  private final BigDecimal expectedMemoryRetryMultiplier = BigDecimal.valueOf(2.0);

  @Test
  void testImputationConfiguration() {
    PipelineConfigurations.ArrayImputationConfig arrayImputationConfiguration =
        pipelineConfigurations.getArrayImputation().get("1");

    assertEquals(1, arrayImputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
    assertEquals(
        "https://test_storage_workspace_url",
        arrayImputationConfiguration.getStorageWorkspaceContainerUrl());
    assertTrue(arrayImputationConfiguration.isUseCallCaching());
    assertFalse(arrayImputationConfiguration.isDeleteIntermediateFiles());
    assertEquals(
        expectedMemoryRetryMultiplier, arrayImputationConfiguration.getMemoryRetryMultiplier());
  }

  @Test
  void testPipelinesCommonConfiguration() {
    PipelineConfigurations.PipelinesCommonConfiguration pipelinesCommonConfiguration =
        pipelineConfigurations.getCommon();
    assertEquals(2, pipelinesCommonConfiguration.getUserDataTtlDays());

    assertEquals(1, pipelinesCommonConfiguration.getQuotaConsumedPollingIntervalSeconds());
    assertTrue(pipelinesCommonConfiguration.isQuotaConsumedUseCallCaching());

    assertEquals(1, pipelinesCommonConfiguration.getInputQcPollingIntervalSeconds());
    assertTrue(pipelinesCommonConfiguration.isInputQcUseCallCaching());

    assertEquals(
        "gs://test_bucket/test_path/to/monitoring/script.sh",
        pipelinesCommonConfiguration.getMonitoringScriptPath());
  }
}
