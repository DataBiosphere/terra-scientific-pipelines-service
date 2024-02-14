package bio.terra.pipelines.dependencies.sam;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import org.broadinstitute.dsde.workbench.client.sam.api.StatusApi;
import org.broadinstitute.dsde.workbench.client.sam.api.UsersApi;
import org.broadinstitute.dsde.workbench.client.sam.auth.Authentication;
import org.broadinstitute.dsde.workbench.client.sam.auth.OAuth;
import org.junit.jupiter.api.Test;

class SamClientTest {
  SamClient samClient;

  String expectedBaseUri = "baseuri";
  SamConfiguration samConfiguration = new SamConfiguration(expectedBaseUri);

  @Test
  void testSamAuthorizedClient() {
    samClient = new SamClient(samConfiguration);

    StatusApi statusApi = samClient.statusApi();
    assertEquals(expectedBaseUri, statusApi.getApiClient().getBasePath());

    String expectedToken = "expected_token";
    UsersApi usersApi = samClient.usersApi(expectedToken);
    for (Authentication auth : usersApi.getApiClient().getAuthentications().values()) {
      if (auth instanceof OAuth) {
        String actual_token = ((OAuth) auth).getAccessToken();
        assertEquals(expectedToken, actual_token);
        return; // sam client only adds a token to the first instance in the list
      }
    }
  }
}
