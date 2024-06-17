package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.app.configuration.external.WorkspaceManagerServerConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class WorkspaceManagerServerConfigurationTest extends BaseEmbeddedDbTest {
  /** test reading workspace manager config from application yml */
  @Autowired WorkspaceManagerServerConfiguration workspaceManagerServerConfiguration;

  @Test
  void verifyWorkspaceManagerServerConfig() {
    assertEquals("https://test_workspace_url/", workspaceManagerServerConfiguration.baseUri());
    assertEquals(24L, workspaceManagerServerConfiguration.sasExpirationDurationHours());
    assertTrue(workspaceManagerServerConfiguration.debugApiLogging());
  }
}
