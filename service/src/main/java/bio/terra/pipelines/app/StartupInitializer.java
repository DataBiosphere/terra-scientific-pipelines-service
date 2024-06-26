package bio.terra.pipelines.app;

import bio.terra.common.migrate.LiquibaseMigrator;
import bio.terra.pipelines.app.configuration.external.SentryConfiguration;
import bio.terra.pipelines.app.configuration.internal.TeaspoonsDatabaseConfiguration;
import bio.terra.pipelines.dependencies.stairway.JobService;
import io.sentry.Sentry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.context.ApplicationContext;

public final class StartupInitializer {
  private static final String CHANGELOG_PATH = "db/changelog.xml";
  private static final Logger logger = LoggerFactory.getLogger(StartupInitializer.class);

  public static void initialize(ApplicationContext applicationContext) {
    // Initialize the Terra Scientific Pipelines Service library
    LiquibaseMigrator migrateService = applicationContext.getBean(LiquibaseMigrator.class);
    JobService jobService = applicationContext.getBean(JobService.class);
    TeaspoonsDatabaseConfiguration teaspoonsDatabaseConfiguration =
        applicationContext.getBean(TeaspoonsDatabaseConfiguration.class);

    if (teaspoonsDatabaseConfiguration.isInitializeOnStart()) {
      migrateService.initialize(CHANGELOG_PATH, teaspoonsDatabaseConfiguration.getDataSource());
    } else if (teaspoonsDatabaseConfiguration.isUpgradeOnStart()) {
      migrateService.upgrade(CHANGELOG_PATH, teaspoonsDatabaseConfiguration.getDataSource());
    }

    // The JobService initialization also handles Stairway initialization.
    jobService.initialize();

    // initialize sentry if sentry dsn env variable is available
    SentryConfiguration sentryConfiguration = applicationContext.getBean(SentryConfiguration.class);
    if (sentryConfiguration.dsn().isEmpty()) {
      logger.info("No Sentry DSN found. Starting up without it.");
    } else {
      logger.info("Sentry DSN found. 5xx errors will be sent to Sentry.");
      Sentry.init(
          options -> {
            options.setDsn(sentryConfiguration.dsn());
            options.setEnvironment(sentryConfiguration.environment());
          });
    }
  }
}
