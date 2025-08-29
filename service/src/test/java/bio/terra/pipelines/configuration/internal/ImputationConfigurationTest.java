package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class ImputationConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired private ImputationConfiguration imputationConfiguration;
  private final BigDecimal expectedMemoryRetryMultiplier = BigDecimal.valueOf(2.0);

  @Test
  void testImputationConfiguration() {
    assertEquals(1, imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
    assertEquals(
        List.of("refDict", "referencePanelPathPrefix", "geneticMapsPath"),
        imputationConfiguration.getInputKeysToPrependWithStorageWorkspaceContainerUrl());
    assertEquals(
        "https://test_storage_workspace_url",
        imputationConfiguration.getStorageWorkspaceContainerUrl());
    assertEquals(
        "/test_reference_panel_path_prefix/file_path",
        imputationConfiguration.getInputsWithCustomValues().get("referencePanelPathPrefix"));
    assertTrue(imputationConfiguration.isUseCallCaching());
    assertFalse(imputationConfiguration.isDeleteIntermediateFiles());
    assertFalse(imputationConfiguration.isUseReferenceDisk());
    assertEquals(expectedMemoryRetryMultiplier, imputationConfiguration.getMemoryRetryMultiplier());
  }

  @Test
  void imputationConfigurationWithNullCustomValuesThrows() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl = List.of();
    Map<String, String> inputsWithCustomValuesWithMissingValue =
        Collections.singletonMap("refDict", null);

    assertThrows(
        IllegalArgumentException.class,
        () -> {
          new ImputationConfiguration(
              1L,
              inputKeysToPrependWithStorageWorkspaceContainerUrl,
              "https://test_storage_workspace_url",
              inputsWithCustomValuesWithMissingValue, // this should cause an exception
              true,
              false,
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
          new ImputationConfiguration(
              1L,
              inputKeysToPrependWithStorageWorkspaceContainerUrl,
              "https://test_storage_workspace_url",
              inputsWithCustomValuesWithMissingValue, // this should cause an exception
              true,
              false,
              false,
              expectedMemoryRetryMultiplier);
        });
  }
}
