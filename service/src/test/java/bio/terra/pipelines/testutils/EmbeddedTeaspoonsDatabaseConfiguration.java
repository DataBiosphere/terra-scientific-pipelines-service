package bio.terra.pipelines.testutils;

import bio.terra.pipelines.app.configuration.internal.TeaspoonsDatabaseConfiguration;
import javax.sql.DataSource;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.autoconfigure.condition.ConditionalOnProperty;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Configuration;

@Configuration
@EnableConfigurationProperties
@ConditionalOnProperty(prefix = "datasource", name = "testWithEmbeddedDatabase")
public class EmbeddedTeaspoonsDatabaseConfiguration extends TeaspoonsDatabaseConfiguration {
  @Autowired
  @Qualifier("teaspoonsDataSource")
  private DataSource teaspoonsDataSource;

  @Override
  public DataSource getDataSource() {
    // Lazy allocation of the data source
    return teaspoonsDataSource;
  }
}
