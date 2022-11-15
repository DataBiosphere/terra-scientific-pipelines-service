package bio.terra.pipelines;

import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.ContextConfiguration;

@SpringBootTest
@ActiveProfiles({"test", "human-readable-logging"})
@ContextConfiguration(classes = App.class)
public abstract class BaseSpringBootTest {}
