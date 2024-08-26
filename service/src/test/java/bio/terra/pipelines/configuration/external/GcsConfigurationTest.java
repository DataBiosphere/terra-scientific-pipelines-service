package bio.terra.pipelines.configuration.external;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.GcsConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class GcsConfigurationTest extends BaseEmbeddedDbTest {
  /** test reading gcs config from application yml */
  @Autowired GcsConfiguration gcsConfiguration;

  @Test
  void verifyGcsConfig() {
    assertEquals(24, gcsConfiguration.signedUrlGetDurationHours());
    assertEquals(24, gcsConfiguration.signedUrlPutDurationHours());
  }
}
