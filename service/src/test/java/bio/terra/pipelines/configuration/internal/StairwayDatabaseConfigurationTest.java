package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.StairwayDatabaseConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class StairwayDatabaseConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired StairwayDatabaseConfiguration stairwayDatabaseConfiguration;

  @Test
  void testStairwayDatabaseConfiguration() {
    assertTrue(stairwayDatabaseConfiguration.getForceClean());
    assertFalse(stairwayDatabaseConfiguration.getMigrateUpgrade());
    assertNotNull(stairwayDatabaseConfiguration.getDataSource());
  }
}
