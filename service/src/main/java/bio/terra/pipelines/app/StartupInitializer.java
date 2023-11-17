package bio.terra.pipelines.app;

import bio.terra.common.migrate.LiquibaseMigrator;
import bio.terra.pipelines.app.configuration.internal.TspsDatabaseConfiguration;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.context.ApplicationContext;

public final class StartupInitializer {
  private static final Logger logger = LoggerFactory.getLogger(StartupInitializer.class);
  private static final String CHANGELOG_PATH = "db/changelog.xml";

  public static void initialize(ApplicationContext applicationContext) {
    // Initialize the Terra Scientific Pipelines Service library
    LiquibaseMigrator migrateService = applicationContext.getBean(LiquibaseMigrator.class);
    StairwayJobService stairwayJobService = applicationContext.getBean(StairwayJobService.class);
    TspsDatabaseConfiguration tspsDatabaseConfiguration =
        applicationContext.getBean(TspsDatabaseConfiguration.class);

    if (tspsDatabaseConfiguration.isInitializeOnStart()) {
      migrateService.initialize(CHANGELOG_PATH, tspsDatabaseConfiguration.getDataSource());
    } else if (tspsDatabaseConfiguration.isUpgradeOnStart()) {
      migrateService.upgrade(CHANGELOG_PATH, tspsDatabaseConfiguration.getDataSource());
    }

    // The StairwayJobService initialization also handles Stairway initialization.
    stairwayJobService.initialize();
  }
}
