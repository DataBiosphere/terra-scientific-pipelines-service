package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.util.List;
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
        "https://test_storage_workspace_url/",
        imputationConfiguration.getStorageWorkspaceContainerUrl());
    assertTrue(imputationConfiguration.isUseCallCaching());
    assertFalse(imputationConfiguration.isDeleteIntermediateFiles());
    assertFalse(imputationConfiguration.isUseReferenceDisk());
  }
}
