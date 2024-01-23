package bio.terra.pipelines.testutils;

import io.zonky.test.db.provider.postgres.PostgreSQLContainerCustomizer;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

@Configuration
public class EmbeddedPostgresConfiguration {

  @Bean
  public PostgreSQLContainerCustomizer postgresContainerCustomizer() {
    return container -> container.withInitScript("init-embedded-postgres-db.sql");
  }
}
