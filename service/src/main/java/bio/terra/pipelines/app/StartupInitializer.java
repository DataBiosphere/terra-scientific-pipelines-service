package bio.terra.pipelines.app;

import bio.terra.common.migrate.LiquibaseMigrator;
import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import org.springframework.context.ApplicationContext;

public final class StartupInitializer {
  private static final String CHANGELOG_PATH = "db/changelog.xml";

  public static void initialize(ApplicationContext applicationContext) {
    // Initialize the Terra Scientific Pipelines Service library
    LiquibaseMigrator migrateService = applicationContext.getBean(LiquibaseMigrator.class);

    TspsDatabaseConfiguration tspsDatabaseConfiguration =
        applicationContext.getBean(TspsDatabaseConfiguration.class);

    if (tspsDatabaseConfiguration.isInitializeOnStart()) {
      migrateService.initialize(CHANGELOG_PATH, tspsDatabaseConfiguration.getDataSource());
    } else if (tspsDatabaseConfiguration.isUpgradeOnStart()) {
      migrateService.upgrade(CHANGELOG_PATH, tspsDatabaseConfiguration.getDataSource());
    }
  }
}
