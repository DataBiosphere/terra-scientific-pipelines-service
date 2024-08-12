package bio.terra.pipelines;

import bio.terra.common.logging.LoggingInitializer;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.boot.context.properties.ConfigurationPropertiesScan;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.orm.jpa.JpaTransactionManager;
import org.springframework.retry.annotation.EnableRetry;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;
import org.thymeleaf.templatemode.TemplateMode;
import org.thymeleaf.templateresolver.ClassLoaderTemplateResolver;

@SpringBootApplication()
@ComponentScan(
    basePackages = {
      // Dependencies for Stairway
      "bio.terra.common.kubernetes",
      // Scan for iam components & configs
      "bio.terra.common.iam",
      // Scan for logging-related components & configs
      "bio.terra.common.logging",
      // Scan for Liquibase migration components & configs
      "bio.terra.common.migrate",
      // Transaction management and DB retry configuration
      "bio.terra.common.retry.transaction",
      // Stairway initialization and status
      "bio.terra.common.stairway",
      // Scan all Teaspoons service packages
      "bio.terra.pipelines",
    })
@ConfigurationPropertiesScan("bio.terra.pipelines")
@EnableRetry
@EnableTransactionManagement
public class App {
  public static void main(String[] args) {
    System.setProperty("spring.config.name", "application");

    new SpringApplicationBuilder(App.class).initializers(new LoggingInitializer()).run(args);
  }

  // This bean plus the @EnableTransactionManagement annotation above enables the use of the
  // @Transaction annotation to control the transaction properties of the data source.
  @Bean("transactionManager")
  public PlatformTransactionManager getTransactionManager() {
    return new JpaTransactionManager();
  }

  /**
   * This bean is used to resolve the location of the Thymeleaf templates that are used to generate
   * the OpenAPI documentation. The default resolver is used to resolve the location of the Swagger
   * UI index.html file in templates/.
   */
  public ClassLoaderTemplateResolver secondaryTemplateResolver() {
    ClassLoaderTemplateResolver secondaryTemplateResolver = new ClassLoaderTemplateResolver();
    secondaryTemplateResolver.setPrefix("static/");
    secondaryTemplateResolver.setSuffix(".yml");
    secondaryTemplateResolver.setTemplateMode(TemplateMode.TEXT);
    secondaryTemplateResolver.setCharacterEncoding("UTF-8");
    secondaryTemplateResolver.setOrder(1);
    secondaryTemplateResolver.setCheckExistence(true);

    return secondaryTemplateResolver;
  }
}
