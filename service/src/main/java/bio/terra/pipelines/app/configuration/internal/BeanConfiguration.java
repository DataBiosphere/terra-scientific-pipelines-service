package bio.terra.pipelines.app.configuration.internal;

import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.BearerTokenFactory;
import jakarta.servlet.http.HttpServletRequest;
import javax.sql.DataSource;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.jdbc.support.JdbcTransactionManager;
import org.springframework.retry.annotation.EnableRetry;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;
import org.springframework.web.context.annotation.RequestScope;

@Configuration
@EnableRetry
@EnableTransactionManagement
@EnableConfigurationProperties
public class BeanConfiguration {

  private final DataSource dataSource;

  public BeanConfiguration(DataSource dataSource) {
    this.dataSource = dataSource;
  }

  /**
   * Taken from <a
   * href="https://github.com/DataBiosphere/terra-data-catalog/blob/5cda83aef8548ff98e7cfbe2a6eaaed9ad1bff45/common/src/main/java/bio/terra/catalog/config/BeanConfig.java#L34-L38">Terra
   * Data Catalog</a> Lasts for the duration of one request, and is injected into dependent beans,
   * even singletons
   */
  @Bean("bearerToken")
  @RequestScope
  public BearerToken bearerToken(HttpServletRequest request) {
    return new BearerTokenFactory().from(request);
  }

  // This bean plus the @EnableTransactionManagement annotation above enables the use of the
  // @Transaction annotation to control the transaction properties of the data source.
  @Bean("transactionManager")
  public PlatformTransactionManager getTransactionManager() {
    return new JdbcTransactionManager(this.dataSource);
  }
}
