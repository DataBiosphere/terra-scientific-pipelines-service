package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.internal.OidcConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class OidcConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired private OidcConfiguration oidcConfiguration;

  @Test
  void testOidcConfiguration() {
    assertEquals("test_client_id", oidcConfiguration.clientId());
    assertEquals("test_authority_endpoint", oidcConfiguration.authorityEndpoint());
    assertEquals("test_token_endpoint", oidcConfiguration.tokenEndpoint());
  }
}
