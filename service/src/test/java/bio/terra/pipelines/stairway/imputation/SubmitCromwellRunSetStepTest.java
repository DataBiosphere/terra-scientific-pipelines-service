package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.MethodDetails;
import bio.terra.cbas.model.MethodListResponse;
import bio.terra.cbas.model.MethodVersionDetails;
import bio.terra.cbas.model.RunSetRequest;
import bio.terra.cbas.model.RunSetStateResponse;
import bio.terra.cbas.model.WorkflowInputDefinition;
import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
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
  }

  @Test
  void doStepSuccess() throws InterruptedException {
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
    when(cbasService.getAllMethods(any(), any())).thenReturn(getAllMethodsResponse);
    when(cbasService.createRunSet(any(), any(), runSetRequestCaptor.capture()))
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
    when(pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(any(), any()))
        .thenReturn(testWorkflowInputDefinitions);

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // extract the captured RunSetRequest
    RunSetRequest capturedRunSetRequest = runSetRequestCaptor.getValue();
    assertFalse(capturedRunSetRequest.isCallCachingEnabled());
    List<WorkflowInputDefinition> capturedWorkflowInputDefinitions =
        capturedRunSetRequest.getWorkflowInputDefinitions();
    // make sure the workflow definitions are populated
    assertEquals(testWorkflowInputDefinitions.size(), capturedWorkflowInputDefinitions.size());

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        runSetId,
        flightContext.getWorkingMap().get(RunImputationJobFlightMapKeys.RUN_SET_ID, UUID.class));
  }

  @Test
  void doStepNoMatchingMethod() throws InterruptedException {
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
    when(cbasService.getAllMethods(any(), any())).thenReturn(getAllMethodsResponse);

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // make sure the step was a fatal faiilure
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
    when(cbasService.getAllMethods(any(), any())).thenReturn(getAllMethodsResponse);
    when(cbasService.createRunSet(any(), any(), any()))
        .thenThrow(new CbasServiceApiException("cbas error"));

    // do the step, expect a RetryException
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    assertThrows(RetryException.class, () -> submitCromwellRunSetStep.doStep(flightContext));
  }

  @Test
  void undoStepSuccess() throws InterruptedException {
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService, pipelinesService, cbasConfiguration);
    StepResult result = submitCromwellRunSetStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
