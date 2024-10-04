package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestFlightContext;
import bio.terra.stairway.FlightStatus;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class StairwayFailedMetricsCounterHookTest extends BaseEmbeddedDbTest {
  StairwayFailedMetricsCounterHook stairwayFailedMetricsCounterHook =
      new StairwayFailedMetricsCounterHook();

  private SimpleMeterRegistry meterRegistry;

  private final String meterRegistryName = "teaspoons.pipeline.failed.count";

  @BeforeEach
  void setup() {
    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.add(meterRegistry);
    meterRegistry.clear();
  }

  @Test
  void endFlight_notPipelineRunTypeFlight_success() throws InterruptedException {
    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.testing.flight.TestFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter

    stairwayFailedMetricsCounterHook.startFlight(context);
    stairwayFailedMetricsCounterHook.endFlight(context);

    // logic should not be executed because this TestFlight does not contain the
    // do_set_pipeline_run_status_failed_hook key in the inputParameters
    Counter counter = meterRegistry.find(meterRegistryName).counter();
    assertNull(counter);
  }

  @Test
  void endFlight_pipelineRunTypeFlight_success() throws InterruptedException {

    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.pipelines.stairway.imputation.RunImputationGcpJobFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter

    // this includes setting the DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK key to true
    StairwayTestUtils.constructCreateJobInputs(context.getInputParameters());

    stairwayFailedMetricsCounterHook.startFlight(context);

    context.flightStatus(FlightStatus.SUCCESS);

    stairwayFailedMetricsCounterHook.endFlight(context);

    // the flight did not fail, so the metric should not be incremented
    Counter counter = meterRegistry.find(meterRegistryName).counter();
    assertNull(counter);
  }

  @Test
  void endFlight_notPipelineRunTypeFlight_error() throws InterruptedException {
    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.testing.flight.TestFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter

    stairwayFailedMetricsCounterHook.startFlight(context);
    // make the flight fail
    context.flightStatus(FlightStatus.ERROR);
    stairwayFailedMetricsCounterHook.endFlight(context);

    // logic should not be executed because this TestFlight does not contain the
    // do_set_pipeline_run_status_failed_hook key in the inputParameters
    Counter counter = meterRegistry.find(meterRegistryName).counter();
    assertNull(counter);
  }

  @Test
  void endFlight_PipelineRunTypeFlight_error() throws InterruptedException {

    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.pipelines.stairway.imputation.RunImputationGcpJobFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter

    // this includes setting the DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK key to true
    StairwayTestUtils.constructCreateJobInputs(context.getInputParameters());

    stairwayFailedMetricsCounterHook.startFlight(context);

    // make the flight fail
    context.flightStatus(FlightStatus.ERROR);
    assertEquals(FlightStatus.ERROR, context.getFlightStatus());

    stairwayFailedMetricsCounterHook.endFlight(context);

    // should have incremented the metric
    Counter counter = meterRegistry.find(meterRegistryName).counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }
}
