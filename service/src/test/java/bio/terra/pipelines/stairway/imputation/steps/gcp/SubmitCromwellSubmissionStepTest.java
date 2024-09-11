package bio.terra.pipelines.stairway.imputation.steps.gcp;

import static bio.terra.pipelines.testutils.TestUtils.VALID_METHOD_CONFIGURATION;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.*;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Captor;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class SubmitCromwellSubmissionStepTest extends BaseEmbeddedDbTest {
  @Mock private RawlsService rawlsService;
  @Captor private ArgumentCaptor<SubmissionRequest> submissionRequestCaptor;
  @Captor private ArgumentCaptor<MethodConfiguration> updateMethodConfigCaptor;
  @Captor private ArgumentCaptor<MethodConfiguration> setMethodConfigCaptor;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;
  @Autowired private ImputationConfiguration imputationConfiguration;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.getCurrentMethodConfigForMethod(any(), any(), any(), any()))
        .thenReturn(new MethodConfiguration());
    when(rawlsService.validateMethodConfig(any(), any(), any(), any(), any(), any()))
        .thenReturn(true);
    when(rawlsService.submitWorkflow(any(), submissionRequestCaptor.capture(), any(), any()))
        .thenReturn(new SubmissionReport().submissionId(testJobId.toString()));

    // do the step
    SubmitCromwellSubmissionStep submitCromwellSubmissionStep =
        new SubmitCromwellSubmissionStep(rawlsService, samService, imputationConfiguration);
    StepResult result = submitCromwellSubmissionStep.doStep(flightContext);

    // extract the captured RunSetRequest and validate
    SubmissionRequest submissionRequest = submissionRequestCaptor.getValue();
    assertFalse(submissionRequest.isDeleteIntermediateOutputFiles());
    assertFalse(submissionRequest.isUseReferenceDisks());
    assertTrue(submissionRequest.isUseCallCache());
    assertEquals(testJobId.toString(), submissionRequest.getEntityName());

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        testJobId,
        flightContext.getWorkingMap().get(RunImputationJobFlightMapKeys.SUBMISSION_ID, UUID.class));
  }

  @Test
  void doStepWithInvalidMethodConfig() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    // set up "current" method config to be version 1.1.1 with the corresponding method uri
    when(rawlsService.getCurrentMethodConfigForMethod(any(), any(), any(), any()))
        .thenReturn(
            new MethodConfiguration()
                .methodRepoMethod(
                    new MethodRepoMethod()
                        .methodUri("http/path/to/wdl/1.1.1")
                        .methodVersion("1.1.1")));
    when(rawlsService.validateMethodConfig(any(), any(), any(), any(), any(), any()))
        .thenReturn(false);
    when(rawlsService.updateMethodConfigToBeValid(
            updateMethodConfigCaptor.capture(), any(), any(), any(), any(), any()))
        .thenReturn(VALID_METHOD_CONFIGURATION);
    when(rawlsService.setMethodConfigForMethod(
            any(), setMethodConfigCaptor.capture(), any(), any(), any()))
        .thenReturn(null);
    when(rawlsService.submitWorkflow(any(), any(), any(), any()))
        .thenReturn(new SubmissionReport().submissionId(testJobId.toString()));

    // do the step
    SubmitCromwellSubmissionStep submitCromwellSubmissionStep =
        new SubmitCromwellSubmissionStep(rawlsService, samService, imputationConfiguration);
    submitCromwellSubmissionStep.doStep(flightContext);

    // extract the captured updateMethodConfig input and setMethodConfig input and validate
    MethodConfiguration updatedMethodConfigInput = updateMethodConfigCaptor.getValue();
    assertEquals("1.1.1", updatedMethodConfigInput.getMethodRepoMethod().getMethodVersion());
    assertEquals(
        "http/path/to/wdl/1.1.1", updatedMethodConfigInput.getMethodRepoMethod().getMethodUri());

    MethodConfiguration setMethodConfigInput = setMethodConfigCaptor.getValue();
    assertEquals(VALID_METHOD_CONFIGURATION, setMethodConfigInput);
  }

  @Test
  void doStepRawlsErrorRetry() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.getCurrentMethodConfigForMethod(any(), any(), any(), any()))
        .thenReturn(new MethodConfiguration());
    when(rawlsService.validateMethodConfig(any(), any(), any(), any(), any(), any()))
        .thenReturn(true);

    // throw exception on submitting workflow
    when(rawlsService.submitWorkflow(any(), any(), any(), any()))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));
    // do the step
    SubmitCromwellSubmissionStep submitCromwellSubmissionStep =
        new SubmitCromwellSubmissionStep(rawlsService, samService, imputationConfiguration);
    StepResult result = submitCromwellSubmissionStep.doStep(flightContext);
    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());

    // throw exception on setting method config to be valid
    when(rawlsService.validateMethodConfig(any(), any(), any(), any(), any(), any()))
        .thenReturn(false);
    when(rawlsService.setMethodConfigForMethod(any(), any(), any(), any(), any()))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));
    // do the step
    submitCromwellSubmissionStep =
        new SubmitCromwellSubmissionStep(rawlsService, samService, imputationConfiguration);
    result = submitCromwellSubmissionStep.doStep(flightContext);
    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());

    // throw exception getting method config
    when(rawlsService.getCurrentMethodConfigForMethod(any(), any(), any(), any()))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));
    // do the step
    submitCromwellSubmissionStep =
        new SubmitCromwellSubmissionStep(rawlsService, samService, imputationConfiguration);
    result = submitCromwellSubmissionStep.doStep(flightContext);
    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    SubmitCromwellSubmissionStep submitCromwellSubmissionStep =
        new SubmitCromwellSubmissionStep(rawlsService, samService, imputationConfiguration);
    StepResult result = submitCromwellSubmissionStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
