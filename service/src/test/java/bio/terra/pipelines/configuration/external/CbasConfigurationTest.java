package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class CbasConfigurationTest extends BaseEmbeddedDbTest {
  /** test reading Cbas config from application yml */
  @Autowired CbasConfiguration cbasConfiguration;

  @Test
  void verifyCbasServerConfig() {
    assertTrue(cbasConfiguration.debugApiLogging());
  }
}
