package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.time.Duration;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppType;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class LeonardoServerConfigurationTest extends BaseContainerTest {
  @Autowired LeonardoServerConfiguration leonardoServerConfiguration;

  /** test reading leonardo server config from application yml */
  @Test
  void verifyLeonardoServerConfig() {
    assertEquals("https://test_leonardo_url/", leonardoServerConfiguration.baseUri());
    assertEquals(List.of(AppType.CROMWELL), leonardoServerConfiguration.cbasAppTypeNames());
    assertEquals(List.of(AppType.WDS), leonardoServerConfiguration.wdsAppTypeNames());
    assertEquals(Duration.ofSeconds(300), leonardoServerConfiguration.dependencyUrlCacheTtl());
    assertFalse(leonardoServerConfiguration.debugApiLogging());
  }
}
