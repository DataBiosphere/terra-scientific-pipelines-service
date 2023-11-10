package bio.terra.pipelines.configuration;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.time.Duration;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppType;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class LeonardoConfigurationTest extends BaseContainerTest {
  @Autowired LeonardoServerConfiguration leonardoServerConfiguration;

  /** test reading leonardo server config from application yml */
  @Test
  void verifyLeonardoServerConfig() {
    assertEquals(leonardoServerConfiguration.baseUri(), "https://test_leonardo_url/");
    assertEquals(leonardoServerConfiguration.cromwellAppTypeNames(), List.of(AppType.CROMWELL));
    assertEquals(leonardoServerConfiguration.wdsAppTypeNames(), List.of(AppType.WDS));
    assertEquals(leonardoServerConfiguration.dependencyUrlCacheTtl(), Duration.ofSeconds(300));
    assertFalse(leonardoServerConfiguration.debugApiLogging());
  }
}
