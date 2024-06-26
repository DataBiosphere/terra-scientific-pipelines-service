package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.assertFalse;

import bio.terra.pipelines.app.configuration.internal.TeaspoonsDatabaseConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class TeaspoonsDatabaseConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired TeaspoonsDatabaseConfiguration teaspoonsDatabaseConfiguration;

  @Test
  void verifyTeaspoonsDatabaseConfiguration() {
    teaspoonsDatabaseConfiguration.setInitializeOnStart(false);
    teaspoonsDatabaseConfiguration.setUpgradeOnStart(false);
    assertFalse(teaspoonsDatabaseConfiguration.isInitializeOnStart());
    assertFalse(teaspoonsDatabaseConfiguration.isUpgradeOnStart());
  }
}
