package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.List;
import java.util.Set;
import java.util.UUID;
import java.util.stream.Collectors;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class RunJobFlightTest extends BaseEmbeddedDbTest {

  @Autowired private JobService jobService;

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final PipelinesEnum imputationPipelineName = PipelinesEnum.IMPUTATION_BEAGLE;

  private static final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private static final String testUserId = TestUtils.TEST_USER_ID_1;

  private static final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;

  private final List<String> expectedStepNames =
      List.of(
          "AddWdsRowStep",
          "CheckCbasHealthStep",
          "CheckLeonardoHealthStep",
          "CheckWdsHealthStep",
          "GetAppUrisStep",
          "PollCromwellRunSetStatusStep",
          "SubmitCromwellRunSetStep",
          "WriteJobToDbStep");

  @Autowired FlightBeanBag flightBeanBag;
  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void setup() {
    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.add(meterRegistry);
  }

  @AfterEach
  void tearDown() {
    meterRegistry.clear();
    Metrics.globalRegistry.clear();
  }

  @Test
  void createJobFlightSetup() {
    // this tests the setters for this flight in JobBuilder
    assertDoesNotThrow(
        () ->
            jobService
                .newJob()
                .jobId(testJobId)
                .flightClass(RunImputationJobFlight.class)
                .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), "test RunImputationJobFlight")
                .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
                .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName)
                .addParameter(RunImputationJobFlightMapKeys.PIPELINE_ID, testPipelineId)
                .addParameter(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, testPipelineInputs));
  }

  @Test
  void expectedStepsInFlight() {
    RunImputationJobFlight runImputationJobFlight =
        new RunImputationJobFlight(StairwayTestUtils.CREATE_JOB_INPUT_PARAMS, flightBeanBag);
    assertEquals(expectedStepNames.size(), runImputationJobFlight.getSteps().size());

    Set<String> stepNames =
        runImputationJobFlight.getSteps().stream()
            .map(step -> step.getClass().getSimpleName())
            .collect(Collectors.toSet());
    for (String step : expectedStepNames) {
      assertTrue(stepNames.contains(step));
    }

    Counter counter = meterRegistry.find("tsps.pipeline.run.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }

  @Test
  void pipelineRunCountIncremented() {
    Counter counter = meterRegistry.find("tsps.pipeline.run.count").counter();
    assertNull(counter);

    // run setup so counter gets incremented
    new RunImputationJobFlight(StairwayTestUtils.CREATE_JOB_INPUT_PARAMS, flightBeanBag);

    counter = meterRegistry.find("tsps.pipeline.run.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }
}
