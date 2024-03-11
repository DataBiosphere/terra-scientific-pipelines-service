package bio.terra.pipelines.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.testutils.BaseTest;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class MetricsUtilsTest extends BaseTest {

  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void setUp() {
    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.add(meterRegistry);
  }

  @AfterEach
  void tearDown() {
    meterRegistry.clear();
    Metrics.globalRegistry.clear();
  }

  @Test
  void incrementImputationPipelineRunMetrics() {
    PipelinesEnum pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4;

    // increment counter once
    MetricsUtils.incrementPipelineRun(pipelineName);

    Counter counter = meterRegistry.find("tsps.pipeline.run.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
    assertEquals(pipelineName.getValue(), counter.getId().getTag("pipeline"));

    // increment counter again
    MetricsUtils.incrementPipelineRun(pipelineName);
    counter = meterRegistry.find("tsps.pipeline.run.count").counter();
    assertEquals(2, counter.count());
  }

  @Test
  void incrementInputationPipelineFailedMetrics() {
    PipelinesEnum pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4;

    // increment counter once
    MetricsUtils.incrementPipelineRunFailed(pipelineName);

    Counter counter = meterRegistry.find("tsps.pipeline.failed.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
    assertEquals(pipelineName.getValue(), counter.getId().getTag("pipeline"));

    // increment counter again
    MetricsUtils.incrementPipelineRunFailed(pipelineName);
    counter = meterRegistry.find("tsps.pipeline.failed.count").counter();
    assertEquals(2, counter.count());
  }
}
