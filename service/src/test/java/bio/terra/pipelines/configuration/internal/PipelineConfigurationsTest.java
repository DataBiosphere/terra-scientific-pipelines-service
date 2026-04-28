package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelineConfigurationsTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineConfigurations pipelineConfigurations;
  private final BigDecimal expectedMemoryRetryMultiplier = BigDecimal.valueOf(2.0);

  @Test
  void testImputationConfiguration() {
    PipelineConfigurations.WdlBasedPipelineConfig wdlBasedPipelineConfiguration =
        pipelineConfigurations.getArrayImputation().get("1");

    assertEquals(1, wdlBasedPipelineConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
    assertEquals(
        List.of("refDict", "referencePanelPathPrefix", "geneticMapsPath"),
        wdlBasedPipelineConfiguration.getInputKeysToPrependWithStorageWorkspaceContainerUrl());
    assertEquals(
        "https://test_storage_workspace_url",
        wdlBasedPipelineConfiguration.getStorageWorkspaceContainerUrl());
    assertEquals(
        "/test_reference_panel_path_prefix/file_path",
        wdlBasedPipelineConfiguration.getInputsWithCustomValues().get("referencePanelPathPrefix"));
    assertTrue(wdlBasedPipelineConfiguration.isUseCallCaching());
    assertFalse(wdlBasedPipelineConfiguration.isDeleteIntermediateFiles());
    assertEquals(
        expectedMemoryRetryMultiplier, wdlBasedPipelineConfiguration.getMemoryRetryMultiplier());
  }

  @Test
  void imputationConfigurationWithNullCustomValuesThrows() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl = List.of();
    Map<String, String> inputsWithCustomValuesWithMissingValue =
        Collections.singletonMap("refDict", null);

    assertThrows(
        IllegalArgumentException.class,
        () -> {
          new PipelineConfigurations.WdlBasedPipelineConfig(
              1L,
              inputKeysToPrependWithStorageWorkspaceContainerUrl,
              "https://test_storage_workspace_url",
              inputsWithCustomValuesWithMissingValue, // this should cause an exception
              true,
              false,
              expectedMemoryRetryMultiplier);
        });
  }

  @Test
  void imputationConfigurationWithBlankCustomValuesThrows() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl = List.of();
    Map<String, String> inputsWithCustomValuesWithMissingValue =
        Collections.singletonMap("refDict", "");

    assertThrows(
        IllegalArgumentException.class,
        () -> {
          new PipelineConfigurations.WdlBasedPipelineConfig(
              1L,
              inputKeysToPrependWithStorageWorkspaceContainerUrl,
              "https://test_storage_workspace_url",
              inputsWithCustomValuesWithMissingValue, // this should cause an exception
              true,
              false,
              expectedMemoryRetryMultiplier);
        });
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
