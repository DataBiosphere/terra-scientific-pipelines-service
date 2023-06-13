package bio.terra.pipelines.db;

import bio.terra.pipelines.app.StartupInitializer;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.TestInstance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.orm.jpa.DataJpaTest;
import org.springframework.context.ApplicationContext;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.transaction.annotation.Transactional;

@DataJpaTest(properties = "spring.main.lazy-initialization=true")
@Transactional
@ActiveProfiles({"test", "human-readable-logging"})
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
public abstract class BaseDaoTest {
  @Autowired ApplicationContext applicationContext;

  // The initialize call kicks off liquibase to rebuild the db.  Ideally, it's done once at the
  // beginning of each set of tests with the @BeforeAll annotation.  But due to transactions not
  // actually rolling back the data between tests, I have given it a @BeforeEach so it rebuilds
  // the db between every test.  This is cheap because it's an in-memory db, so the cost isn't
  // too high.
  // The @TestInstance(TestInstance.Lifecycle.PER_CLASS) annotation on the class above is to
  // allow the @BeforeAll annotation to be applied to a non-static class, which is necessary
  // to use the injected ApplicationContext
  //  @BeforeAll
  @BeforeEach
  void setupDatabase() {
    StartupInitializer.initialize(applicationContext);
  }
}
