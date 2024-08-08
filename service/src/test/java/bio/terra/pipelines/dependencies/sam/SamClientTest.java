package bio.terra.pipelines.dependencies.sam;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import io.opentelemetry.api.OpenTelemetry;
import org.broadinstitute.dsde.workbench.client.sam.api.StatusApi;
import org.broadinstitute.dsde.workbench.client.sam.api.UsersApi;
import org.broadinstitute.dsde.workbench.client.sam.auth.Authentication;
import org.broadinstitute.dsde.workbench.client.sam.auth.OAuth;
import org.junit.jupiter.api.Test;

class SamClientTest {
  SamClient samClient;

  String expectedBaseUri = "baseuri";
  SamConfiguration samConfiguration = new SamConfiguration(expectedBaseUri);
  OpenTelemetry openTelemetry = OpenTelemetry.noop();

  @Test
  void testSamAuthorizedClient() {
    samClient = new SamClient(samConfiguration, openTelemetry);

    StatusApi statusApi = samClient.statusApi();
    assertEquals(expectedBaseUri, statusApi.getApiClient().getBasePath());

    String expectedToken = "expected_token";
    UsersApi usersApi = samClient.usersApi(expectedToken);
    for (Authentication auth : usersApi.getApiClient().getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actualToken = ((OAuth) auth).getAccessToken();
        assertEquals(expectedToken, actualToken);
        return; // sam client only adds a token to the first instance in the list
      }
    }
  }
}
