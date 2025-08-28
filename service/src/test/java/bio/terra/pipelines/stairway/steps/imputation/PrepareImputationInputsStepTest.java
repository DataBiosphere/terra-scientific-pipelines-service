package bio.terra.pipelines.stairway.steps.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class PrepareImputationInputsStepTest extends BaseEmbeddedDbTest {

  @Mock PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired PipelinesService pipelinesService;
  @Autowired PipelinesRepository pipelinesRepository;
  @Autowired ImputationConfiguration imputationConfiguration;
  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.add(meterRegistry);
  }

  @AfterEach
  void tearDown() {
    meterRegistry.clear();
    Metrics.globalRegistry.clear();
  }

  @Test
  void doStepSuccess() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(
        flightContext.getInputParameters(),
        PipelinesEnum.ARRAY_IMPUTATION,
        1L,
        TestUtils.TEST_USER_ID_1,
        TestUtils.TEST_PIPELINE_INPUTS_ARRAY_IMPUTATION,
        TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
        TestUtils.CONTROL_WORKSPACE_NAME,
        TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
        TestUtils.GCP_STORAGE_PROTOCOL,
        TestUtils.TEST_DOMAIN,
        TestUtils.TOOL_CONFIG_GENERIC,
        TestUtils.TOOL_CONFIG_GENERIC,
        TestUtils.TOOL_CONFIG_GENERIC);

    // make sure the full inputs are not populated before the step is executed
    assertNull(
        flightContext
            .getWorkingMap()
            .get(ImputationJobMapKeys.ALL_PIPELINE_INPUTS, new TypeReference<>() {}));

    // mock the service call to format the pipeline inputs
    Map<String, Object> expectedFormattedPipelineInputs =
        new HashMap<>(Map.of("doesnt", "matter", "this", "is", "mocked", "data"));
    when(pipelineInputsOutputsService.gatherAndFormatPipelineInputs(
            eq(testJobId),
            eq(TestUtils.TOOL_CONFIG_GENERIC.inputDefinitions()),
            eq(TestUtils.TEST_PIPELINE_INPUTS_ARRAY_IMPUTATION),
            eq(TestUtils.GCP_STORAGE_PROTOCOL + TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME),
            eq(imputationConfiguration.getInputsWithCustomValues()),
            eq(imputationConfiguration.getInputKeysToPrependWithStorageWorkspaceContainerUrl()),
            eq(imputationConfiguration.getStorageWorkspaceContainerUrl())))
        .thenReturn(expectedFormattedPipelineInputs);

    // do the step
    var prepareImputationInputsStep =
        new PrepareImputationInputsStep(pipelineInputsOutputsService, imputationConfiguration);
    var result = prepareImputationInputsStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // get info from the flight context to run checks
    FlightMap workingMap = flightContext.getWorkingMap();

    // make sure the full pipeline inputs were populated in the working map
    Map<String, Object> fullInputs =
        workingMap.get(ImputationJobMapKeys.ALL_PIPELINE_INPUTS, new TypeReference<>() {});
    assertNotNull(fullInputs);
    for (String key : expectedFormattedPipelineInputs.keySet()) {
      assertEquals(expectedFormattedPipelineInputs.get(key), fullInputs.get(key));
    }
  }

  @Test
  void undoStepSuccess() {
    var prepareImputationInputsStep =
        new PrepareImputationInputsStep(pipelineInputsOutputsService, imputationConfiguration);
    var result = prepareImputationInputsStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
