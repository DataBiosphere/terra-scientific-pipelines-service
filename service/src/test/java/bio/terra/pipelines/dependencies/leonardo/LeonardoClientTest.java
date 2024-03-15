package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiClient;
import org.broadinstitute.dsde.workbench.client.leonardo.auth.Authentication;
import org.broadinstitute.dsde.workbench.client.leonardo.auth.OAuth;
import org.junit.jupiter.api.Test;

class LeonardoClientTest {
  LeonardoClient leonardoClient;

  String expectedBaseUri = "baseuri";
  LeonardoServerConfiguration leonardoServerConfiguration =
      new LeonardoServerConfiguration(expectedBaseUri, List.of(), List.of(), 10, true);

  @Test
  void testLeonardoAuthorizedClient() {
    leonardoClient = new LeonardoClient(leonardoServerConfiguration);

    ApiClient apiClient = leonardoClient.getUnauthorizedApiClient();

    assertEquals(expectedBaseUri, apiClient.getBasePath());
    assertTrue(apiClient.isDebugging());

    String expectedToken = "expected_token";
    apiClient = leonardoClient.getApiClient(expectedToken);
    for (Authentication auth : apiClient.getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actual_token = ((OAuth) auth).getAccessToken();
        assertEquals(expectedToken, actual_token);
        return; // leonardo client only adds a token to the first instance in the list
      }
    }
  }
}
