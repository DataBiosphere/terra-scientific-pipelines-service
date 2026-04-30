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

  @Test
  void testArrayImputationV1Configuration() {
    // note these are the values in pipelines-config-test.yml and not production values
    PipelineConfigurations.WdlBasedPipelineConfig wdlBasedPipelineConfiguration =
        pipelineConfigurations.getArrayImputation().get("1");

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(0.0);

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
    assertEquals(memoryRetryMultiplier, wdlBasedPipelineConfiguration.getMemoryRetryMultiplier());
  }

  @Test
  void testArrayImputationV0Configuration() {
    // note these are the values in pipelines-config-test.yml and not production values
    PipelineConfigurations.WdlBasedPipelineConfig wdlBasedPipelineConfiguration =
        pipelineConfigurations.getArrayImputation().get("0");

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(1.4);

    assertEquals(2, wdlBasedPipelineConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
    assertEquals(
        List.of("prepend1", "prepend2"),
        wdlBasedPipelineConfiguration.getInputKeysToPrependWithStorageWorkspaceContainerUrl());
    assertEquals(
        "https://test_storage_workspace_url",
        wdlBasedPipelineConfiguration.getStorageWorkspaceContainerUrl());
    assertEquals(
        "/test_reference_panel_path_prefix/file_path",
        wdlBasedPipelineConfiguration.getInputsWithCustomValues().get("referencePanelPathPrefix"));
    assertFalse(wdlBasedPipelineConfiguration.isUseCallCaching());
    assertFalse(wdlBasedPipelineConfiguration.isDeleteIntermediateFiles());
    assertEquals(memoryRetryMultiplier, wdlBasedPipelineConfiguration.getMemoryRetryMultiplier());
  }

  @Test
  void arrayImputationConfigurationWithNullCustomValuesThrows() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl = List.of();
    Map<String, String> inputsWithCustomValuesWithMissingValue =
        Collections.singletonMap("refDict", null);
    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(0.0);

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
              memoryRetryMultiplier);
        });
  }

  @Test
  void arrayImputationConfigurationWithBlankCustomValuesThrows() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl = List.of();
    Map<String, String> inputsWithCustomValuesWithMissingValue =
        Collections.singletonMap("refDict", "");

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(0.0);

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
              memoryRetryMultiplier);
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
