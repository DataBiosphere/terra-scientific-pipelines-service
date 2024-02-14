package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class IngressConfigurationTest extends BaseEmbeddedDbTest {
  @Autowired IngressConfiguration ingressConfiguration;

  @Test
  void verifyIngressConfiguration() {
    assertEquals(TestUtils.TEST_DOMAIN, ingressConfiguration.getDomainName());
  }
}
