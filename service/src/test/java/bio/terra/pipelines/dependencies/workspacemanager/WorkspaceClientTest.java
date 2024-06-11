package bio.terra.pipelines.dependencies.workspacemanager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.app.configuration.external.WorkspaceServerConfiguration;
import bio.terra.workspace.client.ApiClient;
import bio.terra.workspace.client.auth.Authentication;
import bio.terra.workspace.client.auth.OAuth;
import org.junit.jupiter.api.Test;

class WorkspaceClientTest {
  WorkspaceClient workspaceClient;

  String expectedBaseUri = "baseuri";
  WorkspaceServerConfiguration workspaceServerConfiguration =
      new WorkspaceServerConfiguration(expectedBaseUri, true);

  @Test
  void testWorkspaceAuthorizedClient() {
    workspaceClient = new WorkspaceClient(workspaceServerConfiguration);

    ApiClient apiClient = workspaceClient.getUnauthorizedApiClient();

    assertEquals(expectedBaseUri, apiClient.getBasePath());
    assertTrue(apiClient.isDebugging());

    String expectedToken = "expected_token";
    apiClient = workspaceClient.getApiClient(expectedToken);
    for (Authentication auth : apiClient.getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actual_token = ((OAuth) auth).getAccessToken();
        assertEquals(expectedToken, actual_token);
        return; // workspace client only adds a token to the first instance in the list
      }
    }
  }
}
