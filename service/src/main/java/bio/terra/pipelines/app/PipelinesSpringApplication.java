package bio.terra.pipelines.app;

import bio.terra.common.logging.LoggingInitializer;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.autoconfigure.jdbc.DataSourceAutoConfiguration;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.retry.annotation.EnableRetry;
import org.springframework.transaction.annotation.EnableTransactionManagement;

@SpringBootApplication(
    exclude = {
      // We don't make use of DataSource in this application, so exclude it from scanning.
      DataSourceAutoConfiguration.class,
    })
@ComponentScan(
    basePackages = {
      // Scan for logging-related components & configs
      "bio.terra.common.logging",
      // Scan for Liquibase migration components & configs
      "bio.terra.common.migrate",
      // Transaction management and DB retry configuration
      "bio.terra.common.retry.transaction",
      // Scan all policy service packages
      "bio.terra.pipelines",
    })
@EnableRetry
@EnableTransactionManagement
public class PipelinesSpringApplication {
  public static void main(String[] args) {
    System.setProperty("spring.config.name", "application");

    new SpringApplicationBuilder(PipelinesSpringApplication.class)
        .initializers(new LoggingInitializer())
        .run(args);
  }
}
