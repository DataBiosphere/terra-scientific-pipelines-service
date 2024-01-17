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
  /** Passed to Stairway, true will run the migrate to upgrade the database */
  private boolean migrateUpgrade;
  /**
   * Passed to Stairway, true will drop any existing stairway data and purge the work queue.
   * Otherwise existing flights are recovered.
   */
  private boolean forceClean;

  public boolean getMigrateUpgrade() {
    return migrateUpgrade;
  }

  public void setMigrateUpgrade(boolean migrateUpgrade) {
    this.migrateUpgrade = migrateUpgrade;
  }

  public boolean getForceClean() {
    return forceClean;
  }

  public void setForceClean(boolean forceClean) {
    this.forceClean = forceClean;
  }

  public boolean isMigrateUpgrade() {
    return getMigrateUpgrade();
  }

  public boolean isForceClean() {
    return getForceClean();
  }

  private DataSource dataSource;

  public DataSource getDataSource() {
    // Lazy allocation of the data source
    if (dataSource == null) {
      dataSource = DataSourceInitializer.initializeDataSource(this);
    }
    return dataSource;
  }
}
