package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestFlightContext;
import bio.terra.stairway.FlightStatus;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.stream.Stream;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

class StairwayFailedMetricsCounterHookTest extends BaseEmbeddedDbTest {
  StairwayFailedMetricsCounterHook stairwayFailedMetricsCounterHook =
      new StairwayFailedMetricsCounterHook();

  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void setup() {
    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.clear();
    Metrics.globalRegistry.add(meterRegistry);
  }

  private static Stream<Arguments> flightContexts() {

    return Stream.of(
        // arguments: whether to include hook key in input Params,
        // hook key value, flight status at endFlight time, expected count of the failed metric

        arguments(
            true,
            true,
            FlightStatus.SUCCESS,
            0), // flight was successful, so the metric should not be incremented
        arguments(
            true,
            true,
            FlightStatus.ERROR,
            1), // flight failed, so the metric should be incremented
        arguments(
            true,
            true,
            FlightStatus.FATAL,
            1), // flight failed dismally, so the metric should be incremented
        arguments(
            true,
            true,
            FlightStatus.QUEUED,
            1), // flight in an unexpected status for end of flight, so the metric should be
        // incremented
        arguments(
            false,
            false, // doesn't matter
            FlightStatus.ERROR,
            0), // flight failed, but the hook key was not included in the input parameters
        arguments(
            true, false, FlightStatus.ERROR, 0)); // flight failed, but the hook key value is false
  }

  @ParameterizedTest
  @MethodSource("flightContexts")
  <T> void endFlight(
      boolean includeHookKey,
      boolean hookKeyValue,
      FlightStatus endFlightStatus,
      double expectedMetricCount)
      throws InterruptedException {
    var context = new TestFlightContext();

    if (includeHookKey) {
      // this includes setting the DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK key to true
      StairwayTestUtils.constructCreateJobInputs(context.getInputParameters());
      context
          .getInputParameters()
          .put(JobMapKeys.DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK, hookKeyValue);
    }

    stairwayFailedMetricsCounterHook.startFlight(context);

    // set the end flight status
    context.flightStatus(endFlightStatus);

    stairwayFailedMetricsCounterHook.endFlight(context);

    // should have incremented the metric
    Counter counter = meterRegistry.find("teaspoons.pipeline.failed.count").counter();
    if (expectedMetricCount == 0) {
      assertNull(counter);
    } else {
      assertNotNull(counter);
      assertEquals(expectedMetricCount, counter.count());
    }
  }
}
