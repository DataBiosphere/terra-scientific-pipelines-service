package bio.terra.pipelines.dependencies.rawls;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.api.EntitiesApi;
import bio.terra.rawls.api.StatusApi;
import bio.terra.rawls.api.SubmissionsApi;
import bio.terra.rawls.client.ApiClient;
import bio.terra.rawls.client.auth.Authentication;
import bio.terra.rawls.client.auth.OAuth;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class RawlsClientTest extends BaseEmbeddedDbTest {
  @Autowired RawlsClient rawlsClient;

  String rawlsBaseUri = "https://test_rawls_url/";
  String authToken = "authToken";

  @Test
  void testWorkspaceAuthorizedClient() {
    ApiClient apiClient = rawlsClient.getUnauthorizedApiClient();

    assertEquals(rawlsBaseUri, apiClient.getBasePath());
    assertTrue(apiClient.isDebugging());

    apiClient = rawlsClient.getApiClient(authToken);
    for (Authentication auth : apiClient.getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actualToken = ((OAuth) auth).getAccessToken();
        assertEquals(authToken, actualToken);
        return; // rawls client only adds a token to the first instance in the list
      }
    }
  }

  @Test
  void testWorkspaceClientApis() {
    StatusApi unauthenticatedApi = rawlsClient.getStatusApi();
    assertEquals(rawlsBaseUri, unauthenticatedApi.getApiClient().getBasePath());
    assertTrue(unauthenticatedApi.getApiClient().isDebugging());

    EntitiesApi resourceApi = rawlsClient.getEntitiesApi(authToken);
    assertEquals(rawlsBaseUri, resourceApi.getApiClient().getBasePath());
    assertTrue(resourceApi.getApiClient().isDebugging());

    SubmissionsApi submissionsApi = rawlsClient.getSubmissionsApi(authToken);
    assertEquals(rawlsBaseUri, submissionsApi.getApiClient().getBasePath());
    assertTrue(submissionsApi.getApiClient().isDebugging());
  }
}
