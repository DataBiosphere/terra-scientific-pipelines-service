package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.StairwayDatabaseConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;

@SpringBootTest(
    properties = {
      "spring.main.lazy-initialization=true",
      "datasource.testWithEmbeddedDatabase=false"
    })
class StairwayDatabaseConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired StairwayDatabaseConfiguration stairwayDatabaseConfiguration;

  @Test
  void testStairwayDatabaseConfiguration() {
    assertNotNull(stairwayDatabaseConfiguration.getDataSource());
  }
}
