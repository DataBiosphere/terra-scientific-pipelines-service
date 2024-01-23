package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.assertFalse;

import bio.terra.pipelines.app.configuration.internal.TspsDatabaseConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class TspsDatabaseConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired TspsDatabaseConfiguration tspsDatabaseConfiguration;

  @Test
  void verifyTspsDatabaseConfiguration() {
    assertFalse(tspsDatabaseConfiguration.isInitializeOnStart());
    assertFalse(tspsDatabaseConfiguration.isUpgradeOnStart());
  }
}
