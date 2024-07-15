package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class PrepareImputationInputsStepTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesService pipelinesService;
  @Autowired PipelinesRepository pipelinesRepository;
  @Autowired ImputationConfiguration imputationConfiguration;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);

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

    PipelinesEnum pipelineEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesService.getPipeline(pipelineEnum);

    StairwayTestUtils.constructCreateJobInputs(
        flightContext.getInputParameters(),
        PipelinesEnum.IMPUTATION_BEAGLE,
        pipeline.getId(),
        pipeline.getPipelineInputDefinitions(),
        pipeline.getPipelineOutputDefinitions(),
        TestUtils.TEST_USER_ID_1,
        TestUtils.TEST_PIPELINE_INPUTS_IMPUTATION_BEAGLE,
        TestUtils.CONTROL_WORKSPACE_ID,
        TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_URL,
        pipeline.getWdlMethodName(),
        TestUtils.TEST_RESULT_URL);

    // make sure the full inputs are not populated before the step is executed
    assertNull(
        flightContext
            .getWorkingMap()
            .get(RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, new TypeReference<>() {}));

    // do the step
    var prepareImputationInputsStep =
        new PrepareImputationInputsStep(pipelinesService, imputationConfiguration);
    var result = prepareImputationInputsStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();
    FlightMap workingMap = flightContext.getWorkingMap();

    // get the service-provided inputs
    List<PipelineInputDefinition> serviceProvidedPipelineInputDefinitions =
        pipeline.getPipelineInputDefinitions().stream()
            .filter(Predicate.not(PipelineInputDefinition::getUserProvided))
            .toList();

    // get the user-provided inputs
    List<PipelineInputDefinition> userProvidedPipelineInputDefinitions =
        pipeline.getPipelineInputDefinitions().stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .toList();

    // make sure the full map of inputs was prepared
    Map<String, Object> userProvidedInputs =
        inputParams.get(
            RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS, new TypeReference<>() {});
    Map<String, Object> fullInputs =
        workingMap.get(RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, new TypeReference<>() {});

    // make sure the fullInputs map contains all the user-provided keys as well as all the
    // service-provided keys
    assertNotNull(fullInputs);
    assertNotNull(userProvidedInputs);

    // input definitions use camelCase names, whereas fullInputs use camel_case wdlVariableNames
    for (String wdlInputName :
        userProvidedPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())) {
      assertTrue(fullInputs.containsKey(wdlInputName));
    }
    for (String wdlInputName :
        serviceProvidedPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())) {
      assertTrue(fullInputs.containsKey(wdlInputName));
    }

    // user-provided file inputs should contain the control workspace's storage url
    for (String wdlInputName :
        userProvidedPipelineInputDefinitions.stream()
            .filter(input -> input.getType().equals(PipelineInputTypesEnum.VCF))
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())) {
      assertTrue(
          fullInputs
              .get(wdlInputName)
              .toString()
              .contains(TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_URL));
    }

    // make sure each input in the fullInputs map has a populated value
    for (String inputName : fullInputs.keySet()) {
      assertNotNull(fullInputs.get(inputName));
    }
  }

  @Test
  void undoStepSuccess() {
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    var prepareImputationInputsStep =
        new PrepareImputationInputsStep(pipelinesService, imputationConfiguration);
    var result = prepareImputationInputsStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
