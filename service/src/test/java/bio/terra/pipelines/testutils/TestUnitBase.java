package bio.terra.pipelines.testutils;

import bio.terra.pipelines.app.App;
import org.junit.jupiter.api.extension.ExtendWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit.jupiter.SpringExtension;

@ActiveProfiles({"test", "human-readable-logging"})
@ContextConfiguration(classes = App.class)
@ExtendWith(SpringExtension.class)
@SpringBootTest()
public class TestUnitBase {}
