package bio.terra.pipelines.dependencies.workspacemanager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.workspace.api.ControlledAzureResourceApi;
import bio.terra.workspace.api.ResourceApi;
import bio.terra.workspace.api.UnauthenticatedApi;
import bio.terra.workspace.client.ApiClient;
import bio.terra.workspace.client.auth.Authentication;
import bio.terra.workspace.client.auth.OAuth;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class WorkspaceClientTest extends BaseEmbeddedDbTest {

  @Autowired WorkspaceClient workspaceClient;

  String workspaceBaseUri = "https://test_workspace_url/";
  String authToken = "authToken";

  @Test
  void testWorkspaceAuthorizedClient() {
    ApiClient apiClient = workspaceClient.getUnauthorizedApiClient();

    assertEquals(workspaceBaseUri, apiClient.getBasePath());
    assertTrue(apiClient.isDebugging());

    apiClient = workspaceClient.getApiClient(authToken);
    for (Authentication auth : apiClient.getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actual_token = ((OAuth) auth).getAccessToken();
        assertEquals(authToken, actual_token);
        return; // workspace client only adds a token to the first instance in the list
      }
    }
  }

  @Test
  void testWorkspaceClientApis() {
    UnauthenticatedApi unauthenticatedApi = workspaceClient.getUnauthenticatedApi();
    assertEquals(workspaceBaseUri, unauthenticatedApi.getApiClient().getBasePath());
    assertTrue(unauthenticatedApi.getApiClient().isDebugging());

    ResourceApi resourceApi = workspaceClient.getResourceApi(authToken);
    assertEquals(workspaceBaseUri, resourceApi.getApiClient().getBasePath());
    assertTrue(resourceApi.getApiClient().isDebugging());

    ControlledAzureResourceApi controlledAzureResourceApi =
        workspaceClient.getControlledAzureResourceApi(authToken);
    assertEquals(workspaceBaseUri, controlledAzureResourceApi.getApiClient().getBasePath());
    assertTrue(controlledAzureResourceApi.getApiClient().isDebugging());
  }
}
