package bio.terra.pipelines.testutils;

import bio.terra.pipelines.app.configuration.internal.StairwayDatabaseConfiguration;
import javax.sql.DataSource;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.autoconfigure.condition.ConditionalOnProperty;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Configuration;

@Configuration
@EnableConfigurationProperties
@ConditionalOnProperty(prefix = "datasource", name = "testWithEmbeddedDatabase")
public class EmbeddedStairwayDatabaseConfiguration extends StairwayDatabaseConfiguration {

  @Autowired
  @Qualifier("stairwayDataSource")
  private DataSource stairwayDataSource;

  @Override
  public DataSource getDataSource() {
    // Lazy allocation of the data source
    return stairwayDataSource;
  }
}
