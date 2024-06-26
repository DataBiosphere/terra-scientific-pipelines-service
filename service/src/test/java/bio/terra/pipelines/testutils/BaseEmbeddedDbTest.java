package bio.terra.pipelines.testutils;

import static io.zonky.test.db.AutoConfigureEmbeddedDatabase.RefreshMode.AFTER_EACH_TEST_METHOD;

import io.zonky.test.db.AutoConfigureEmbeddedDatabase;
// import jakarta.activation.DataSource;
import javax.sql.DataSource;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.TestConfiguration;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Primary;
import org.springframework.jdbc.datasource.UserCredentialsDataSourceAdapter;

@SpringBootTest(properties = "spring.main.lazy-initialization=true")
@AutoConfigureEmbeddedDatabase(
    type = AutoConfigureEmbeddedDatabase.DatabaseType.POSTGRES,
    provider = AutoConfigureEmbeddedDatabase.DatabaseProvider.DOCKER,
    beanName = "commonDataSource",
    refresh = AFTER_EACH_TEST_METHOD)
public abstract class BaseEmbeddedDbTest extends BaseTest {

  @TestConfiguration
  static class Config {

    @Bean
    @Primary
    public DataSource teaspoonsDataSource(DataSource commonDataSource) {
      UserCredentialsDataSourceAdapter adapter = new UserCredentialsDataSourceAdapter();
      adapter.setTargetDataSource(commonDataSource);
      adapter.setSchema("public");
      return adapter;
    }

    @Bean
    public DataSource stairwayDataSource(DataSource commonDataSource) {
      UserCredentialsDataSourceAdapter adapter = new UserCredentialsDataSourceAdapter();
      adapter.setTargetDataSource(commonDataSource);
      adapter.setSchema("public");
      return adapter;
    }
  }
}
