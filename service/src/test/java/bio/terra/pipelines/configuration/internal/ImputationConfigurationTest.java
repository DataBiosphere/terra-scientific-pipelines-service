package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class ImputationConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired private ImputationConfiguration imputationConfiguration;

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
  }

  @Test
  void testImputationConfigurationWithInvalidStorageWorkspaceUrl() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl =
        List.of("refDict", "referencePanelPathPrefix", "geneticMapsPath");
    Map<String, String> inputsWithCustomValues = Map.of();
    assertThrows(
        IllegalArgumentException.class,
        () -> {
          new ImputationConfiguration(
              1L,
              inputKeysToPrependWithStorageWorkspaceContainerUrl,
              "https://test_storage_workspace_url/", // this should cause an exception
              inputsWithCustomValues,
              true,
              false,
              false);
        });
  }

  @Test
  void testImputationConfigurationWithMissingInputsWithCustomValues() {
    List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl =
        List.of("refDict", "referencePanelPathPrefix", "geneticMapsPath");

    // need to set this up as a HashMap to allow for null values, which aren't allowed in a Map
    Map<String, String> inputsWithCustomValuesWithMissingValue = new HashMap<>();
    inputsWithCustomValuesWithMissingValue.put("refDict", null);

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
              false);
        });
  }
}
