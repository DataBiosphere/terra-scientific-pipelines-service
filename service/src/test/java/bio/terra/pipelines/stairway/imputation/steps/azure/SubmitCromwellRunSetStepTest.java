package bio.terra.pipelines.stairway.imputation.steps.azure;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.MethodDetails;
import bio.terra.cbas.model.MethodListResponse;
import bio.terra.cbas.model.MethodVersionDetails;
import bio.terra.cbas.model.RunSetRequest;
import bio.terra.cbas.model.RunSetStateResponse;
import bio.terra.cbas.model.WorkflowInputDefinition;
import bio.terra.cbas.model.WorkflowOutputDefinition;
import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.*;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Captor;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class SubmitCromwellRunSetStepTest extends BaseEmbeddedDbTest {
  @Mock private CbasService cbasService;
  @Captor private ArgumentCaptor<RunSetRequest> runSetRequestCaptor;
  @Mock private SamService samService;
  @Mock private PipelinesService pipelinesService;
  @Mock private FlightContext flightContext;
  @Autowired private CbasConfiguration cbasConfiguration;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();
    workingMap.put(RunImputationJobFlightMapKeys.CBAS_URI, "cbasUri");
    workingMap.put(
        RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, TestUtils.TEST_PIPELINE_INPUTS);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    UUID runSetId = UUID.randomUUID();
    MethodListResponse getAllMethodsResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name(
                        flightContext
                            .getInputParameters()
                            .get(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, String.class))
                    .addMethodVersionsItem(
                        new MethodVersionDetails().methodVersionId(UUID.randomUUID())));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.getAllMethods("cbasUri", "thisToken")).thenReturn(getAllMethodsResponse);
    when(cbasService.createRunSet(eq("cbasUri"), eq("thisToken"), runSetRequestCaptor.capture()))
        .thenReturn(new RunSetStateResponse().runSetId(runSetId));

    List<WorkflowInputDefinition> testWorkflowInputDefinitions =
        new ArrayList<>(
            List.of(
                new WorkflowInputDefinition()
                    .inputName("testInputName1")
                    .inputType(null)
                    .source(null),
                new WorkflowInputDefinition()
                    .inputName("testInputName2")
                    .inputType(null)
                    .source(null)));
    when(pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST, TestUtils.TEST_WDL_METHOD_NAME_1))
        .thenReturn(testWorkflowInputDefinitions);

    List<WorkflowOutputDefinition> testWorkflowOutputDefinitions =
        new ArrayList<>(
            List.of(
                new WorkflowOutputDefinition()
                    .outputName("testOutputName1")
                    .outputType(null)
                    .destination(null),
                new WorkflowOutputDefinition()
                    .outputName("testOutputName2")
                    .outputType(null)
                    .destination(null)));
    when(pipelinesService.prepareCbasWorkflowOutputRecordUpdateDefinitions(
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST, TestUtils.TEST_WDL_METHOD_NAME_1))
        .thenReturn(testWorkflowOutputDefinitions);

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // extract the captured RunSetRequest
    RunSetRequest capturedRunSetRequest = runSetRequestCaptor.getValue();
    assertFalse(capturedRunSetRequest.isCallCachingEnabled());
    List<WorkflowInputDefinition> capturedWorkflowInputDefinitions =
        capturedRunSetRequest.getWorkflowInputDefinitions();
    List<WorkflowOutputDefinition> capturedWorkflowOutputDefinitions =
        capturedRunSetRequest.getWorkflowOutputDefinitions();
    // make sure the workflow definitions are populated
    assertEquals(testWorkflowInputDefinitions.size(), capturedWorkflowInputDefinitions.size());
    assertEquals(testWorkflowOutputDefinitions.size(), capturedWorkflowOutputDefinitions.size());

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        runSetId,
        flightContext.getWorkingMap().get(RunImputationJobFlightMapKeys.RUN_SET_ID, UUID.class));
  }

  @Test
  void doStepNoMatchingMethod() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    MethodListResponse getAllMethodsResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name("random name that doesnt match anything")
                    .addMethodVersionsItem(
                        new MethodVersionDetails().methodVersionId(UUID.randomUUID())));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.getAllMethods("cbasUri", "thisToken")).thenReturn(getAllMethodsResponse);

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // make sure the step was a fatal failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }

  @Test
  void doStepCbasErrorRetry() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    MethodListResponse getAllMethodsResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name(
                        flightContext
                            .getInputParameters()
                            .get(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, String.class))
                    .addMethodVersionsItem(
                        new MethodVersionDetails().methodVersionId(UUID.randomUUID())));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.getAllMethods("cbasUri", "thisToken")).thenReturn(getAllMethodsResponse);
    when(cbasService.createRunSet(eq("cbasUri"), eq("thisToken"), any(RunSetRequest.class)))
        .thenThrow(new CbasServiceApiException("cbas error"));

    // do the step, expect a Retry status
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
