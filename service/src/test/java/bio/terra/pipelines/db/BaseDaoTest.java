package bio.terra.pipelines.db;

import bio.terra.pipelines.app.StartupInitializer;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.TestInstance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.data.jdbc.DataJdbcTest;
import org.springframework.context.ApplicationContext;
import org.springframework.test.context.ActiveProfiles;

@DataJdbcTest(properties = "spring.main.lazy-initialization=true")
@ActiveProfiles({"test", "human-readable-logging"})
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
public abstract class BaseDaoTest {
  @Autowired ApplicationContext applicationContext;

  //  @BeforeAll
  @BeforeEach
  void setupDatabase() {
    StartupInitializer.initialize(applicationContext);
  }
}
