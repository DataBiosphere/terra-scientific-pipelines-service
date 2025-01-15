package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelinesCommonConfigurationTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesCommonConfiguration pipelinesCommonConfiguration;

  @Test
  void testPipelinesCommonConfiguration() {
    assertEquals(1, pipelinesCommonConfiguration.quotaConsumedPollingIntervalSeconds());
    assertTrue(pipelinesCommonConfiguration.quotaConsumedUseCallCaching());
    assertEquals(2, pipelinesCommonConfiguration.storageBucketTtlDays());
  }
}
