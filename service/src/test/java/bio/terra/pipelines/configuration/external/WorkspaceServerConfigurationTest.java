package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.app.configuration.external.WorkspaceServerConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class WorkspaceServerConfigurationTest extends BaseEmbeddedDbTest {
  /** test reading workspace manager config from application yml */
  @Autowired WorkspaceServerConfiguration workspaceServerConfiguration;

  @Test
  void verifyWorkspaceServerConfig() {
    assertEquals("https://test_workspace_url/", workspaceServerConfiguration.baseUri());
    assertTrue(workspaceServerConfiguration.debugApiLogging());
  }
}
