package bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.v20260428;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class RunWdlBasedPipelineJobFlightTest extends BaseEmbeddedDbTest {

  @Autowired private JobService jobService;

  private final List<String> expectedStepNames =
      List.of(
          "PrepareInputsStep",
          "AddDataTableRowStep",
          // quota wdl steps
          "SubmitCromwellSubmissionStep",
          "PollCromwellSubmissionStatusStep",
          "FetchOutputsFromDataTableStep",
          "QuotaConsumedValidationStep",
          // input qc wdl steps
          "SubmitCromwellSubmissionStep",
          "PollCromwellSubmissionStatusStep",
          "FetchOutputsFromDataTableStep",
          "InputQcValidationStep",
          // pipeline wdl steps
          "SubmitCromwellSubmissionStep",
          "PollCromwellSubmissionStatusStep",
          "FetchOutputsFromDataTableStep",
          "PopulateFileOutputSizeStep",
          "CompletePipelineRunStep",
          "SendJobSucceededNotificationStep");

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
    // this tests the setters for this flight in JobBuilder. note this doesn't check for required
    // input parameters, nor does it do anything pipeline-specific
    assertDoesNotThrow(
        () ->
            jobService
                .newJob()
                .jobId(TestUtils.TEST_NEW_UUID)
                .flightClass(RunWdlBasedPipelineJobFlight.class)
                .addParameter(JobMapKeys.DESCRIPTION, "test RunWdlBasedPipelineJobFlight")
                .addParameter(JobMapKeys.USER_ID, TestUtils.TEST_USER_1_ID)
                .addParameter(JobMapKeys.PIPELINE_NAME, PipelinesEnum.ARRAY_IMPUTATION)
                .addParameter(JobMapKeys.PIPELINE_VERSION, TestUtils.TEST_PIPELINE_VERSION_1)
                .addParameter(JobMapKeys.PIPELINE_ID, TestUtils.TEST_PIPELINE_ID_1)
                .addParameter(
                    WdlBasedPipelineJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
                    TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST)
                .addParameter(
                    WdlBasedPipelineJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
                    TestUtils.TEST_PIPELINE_INPUTS)
                .addParameter(
                    WdlBasedPipelineJobMapKeys.PIPELINE_TOOL_CONFIG, TestUtils.TOOL_CONFIG_GENERIC)
                .addParameter(
                    WdlBasedPipelineJobMapKeys.QUOTA_TOOL_CONFIG, TestUtils.TOOL_CONFIG_GENERIC)
                .addParameter(
                    WdlBasedPipelineJobMapKeys.INPUT_QC_TOOL_CONFIG,
                    TestUtils.TOOL_CONFIG_GENERIC));
  }

  @Test
  void flightRunsExpectedPipeline() {
    RunWdlBasedPipelineJobFlight runArrayImputationPipelineJobFlight =
        new RunWdlBasedPipelineJobFlight(
            StairwayTestUtils.CREATE_ARRAY_IMPUTATION_JOB_INPUT_PARAMS, flightBeanBag);
    assertEquals(
        PipelinesEnum.ARRAY_IMPUTATION,
        runArrayImputationPipelineJobFlight
            .getInputParameters()
            .get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class));
  }

  @Test
  void expectedStepsInFlight() {
    RunWdlBasedPipelineJobFlight runWdlBasedPipelineJobFlight =
        new RunWdlBasedPipelineJobFlight(
            StairwayTestUtils.CREATE_ARRAY_IMPUTATION_JOB_INPUT_PARAMS, flightBeanBag);
    assertEquals(expectedStepNames.size(), runWdlBasedPipelineJobFlight.getSteps().size());

    Set<String> stepNames =
        runWdlBasedPipelineJobFlight.getSteps().stream()
            .map(step -> step.getClass().getSimpleName())
            .collect(Collectors.toSet());
    for (String step : expectedStepNames) {
      assertTrue(stepNames.contains(step));
    }

    Counter counter = meterRegistry.find("teaspoons.pipeline.run.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }

  @Test
  void pipelineRunCountIncremented() {
    Counter counter = meterRegistry.find("teaspoons.pipeline.run.count").counter();
    assertNull(counter);

    // run setup so counter gets incremented
    new RunWdlBasedPipelineJobFlight(
        StairwayTestUtils.CREATE_ARRAY_IMPUTATION_JOB_INPUT_PARAMS, flightBeanBag);

    counter = meterRegistry.find("teaspoons.pipeline.run.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }
}
