package bio.terra.pipelines.testutils;

import org.testcontainers.containers.PostgreSQLContainer;

public class TestPostgresqlContainer extends PostgreSQLContainer<TestPostgresqlContainer> {
  private static final String IMAGE_VERSION = "postgres:13.1";
  private static TestPostgresqlContainer container;
  private static final String TSPS_DB_NAME = "pipelines_db";
  private static final String STAIRWAY_DB_NAME = "tsps_stairway_db";

  private TestPostgresqlContainer() {
    super(IMAGE_VERSION);
  }

  public static TestPostgresqlContainer getInstance() {
    if (container == null) {
      container =
          new TestPostgresqlContainer()
              .withDatabaseName(TSPS_DB_NAME)
              .withInitScript("init-create-stairway-db.sql");
    }
    return container;
  }

  @Override
  public void start() {
    super.start();
    System.setProperty("TSPS_DB_URL", container.getJdbcUrl());
    System.setProperty("STAIRWAY_DB_URL", getStairwayJdbcUrl());
    System.setProperty("DB_USERNAME", container.getUsername());
    System.setProperty("DB_PASSWORD", container.getPassword());
  }

  @Override
  public void stop() {
    // do nothing, JVM handles shut down
  }

  /*
   * The StairwayJobService uses the same container as the TSPS service, but with a different database name.
   * This method returns the JDBC URL with the correct database name for the Stairway database.
   */
  private String getStairwayJdbcUrl() {
    String fullUrl = container.getJdbcUrl();
    // fullUrl looks like jdbc:postgresql://localhost:50218/pipelines_db?loggerLevel=OFF
    return fullUrl.replace(TSPS_DB_NAME, STAIRWAY_DB_NAME);
  }
}
