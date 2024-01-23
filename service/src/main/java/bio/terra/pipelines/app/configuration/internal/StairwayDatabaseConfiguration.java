package bio.terra.pipelines.app.configuration.internal;

import bio.terra.common.db.BaseDatabaseProperties;
import bio.terra.common.db.DataSourceInitializer;
import javax.sql.DataSource;
import org.springframework.boot.autoconfigure.condition.ConditionalOnProperty;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Configuration;

@Configuration
@EnableConfigurationProperties
@ConfigurationProperties(prefix = "spring.stairway-database")
@ConditionalOnProperty(
    prefix = "datasource",
    name = "testWithEmbeddedDatabase",
    havingValue = "false",
    matchIfMissing = true)
public class StairwayDatabaseConfiguration extends BaseDatabaseProperties {
  private DataSource dataSource;

  public DataSource getDataSource() {
    // Lazy allocation of the data source
    if (dataSource == null) {
      dataSource = DataSourceInitializer.initializeDataSource(this);
    }
    return dataSource;
  }
}
