package bio.terra.pipelines.testutils;

import org.junit.ClassRule;
import org.springframework.boot.test.autoconfigure.jdbc.AutoConfigureTestDatabase;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.transaction.annotation.Transactional;
import org.testcontainers.containers.PostgreSQLContainer;

@Transactional
@SpringBootTest(properties = "spring.main.lazy-initialization=true")
@AutoConfigureTestDatabase(replace = AutoConfigureTestDatabase.Replace.NONE)
public class BaseContainerTest extends BaseTest {
  @ClassRule
  public static PostgreSQLContainer postgreSQLContainer = TestPostgresqlContainer.getInstance();

  static {
    postgreSQLContainer.start();
  }
}
