package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import java.time.Duration;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiClient;
import org.broadinstitute.dsde.workbench.client.leonardo.auth.Authentication;
import org.broadinstitute.dsde.workbench.client.leonardo.auth.OAuth;
import org.junit.jupiter.api.Test;

public class LeonardoClientTest {
  LeonardoClient leonardoClient;
  LeonardoServerConfiguration leonardoServerConfiguration =
      new LeonardoServerConfiguration(
          "baseuri", List.of(), List.of(), Duration.ofMinutes(10), true);

  @Test
  void TestLeonardoUnauthorizedClient() {
    leonardoClient = new LeonardoClient(leonardoServerConfiguration);

    ApiClient apiClient = leonardoClient.getUnauthorizedApiClient();

    assertEquals("baseuri", apiClient.getBasePath());
    assertTrue(apiClient.isDebugging());

    String expectedToken = "expected_token";
    apiClient = leonardoClient.getApiClient(expectedToken);
    for (Authentication auth : apiClient.getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actual_token = ((OAuth) auth).getAccessToken();
        assertEquals(expectedToken, actual_token);
      }
    }
  }
}
