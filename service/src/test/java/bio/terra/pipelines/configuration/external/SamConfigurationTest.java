package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class SamConfigurationTest extends BaseEmbeddedDbTest {
  /** test reading sam config from application yml */
  @Autowired SamConfiguration samConfiguration;

  @Test
  void verifySamServerConfig() {
    assertEquals("testSamUri", samConfiguration.baseUri());
  }
}
