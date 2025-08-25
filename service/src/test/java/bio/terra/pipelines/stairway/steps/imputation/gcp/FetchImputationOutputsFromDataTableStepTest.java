package bio.terra.pipelines.stairway.steps.imputation.gcp;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class FetchImputationOutputsFromDataTableStepTest extends BaseEmbeddedDbTest {

  @Mock RawlsService rawlsService;
  @Mock SamService samService;
  @Mock PipelineInputsOutputsService pipelineInputsOutputsService;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() throws RawlsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    // outputs to match the test output definitions
    Map<String, Object> entityAttributes = new HashMap<>(Map.of("output_name", "some/file.vcf.gz"));
    Entity entity = new Entity().attributes(entityAttributes);

    when(rawlsService.getDataTableEntity(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenReturn(entity);
    Map<String, String> outputsProcessedFromEntity =
        new HashMap<>(Map.of("outputName", "some/file.vcf.gz"));
    when(pipelineInputsOutputsService.extractPipelineOutputsFromEntity(
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST, entity))
        .thenReturn(outputsProcessedFromEntity);

    FetchImputationOutputsFromDataTableStep fetchImputationOutputsFromDataTableStep =
        new FetchImputationOutputsFromDataTableStep(
            rawlsService, samService, pipelineInputsOutputsService);
    StepResult result = fetchImputationOutputsFromDataTableStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    assertEquals(
        outputsProcessedFromEntity,
        flightContext.getWorkingMap().get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class));
  }

  @Test
  void doStepRawlsFailureRetry() throws RawlsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    when(rawlsService.getDataTableEntity(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenThrow(new RawlsServiceApiException("Rawls Service Api Exception"));

    FetchImputationOutputsFromDataTableStep fetchImputationOutputsFromDataTableStep =
        new FetchImputationOutputsFromDataTableStep(
            rawlsService, samService, pipelineInputsOutputsService);
    StepResult result = fetchImputationOutputsFromDataTableStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void doStepOutputsFailureNoRetry() throws InternalServerErrorException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    Map<String, Object> entityAttributes = new HashMap<>(Map.of("output_name", "some/file.vcf.gz"));
    Entity entity = new Entity().attributes(entityAttributes);

    when(rawlsService.getDataTableEntity(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenReturn(entity);
    when(pipelineInputsOutputsService.extractPipelineOutputsFromEntity(
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST, entity))
        .thenThrow(new InternalServerErrorException("Internal Server Error"));

    FetchImputationOutputsFromDataTableStep fetchImputationOutputsFromDataTableStep =
        new FetchImputationOutputsFromDataTableStep(
            rawlsService, samService, pipelineInputsOutputsService);
    assertThrows(
        InternalServerErrorException.class,
        () -> fetchImputationOutputsFromDataTableStep.doStep(flightContext));
  }

  @Test
  void undoStepSuccess() {
    FetchImputationOutputsFromDataTableStep fetchImputationOutputsFromDataTableStep =
        new FetchImputationOutputsFromDataTableStep(
            rawlsService, samService, pipelineInputsOutputsService);
    StepResult result = fetchImputationOutputsFromDataTableStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
