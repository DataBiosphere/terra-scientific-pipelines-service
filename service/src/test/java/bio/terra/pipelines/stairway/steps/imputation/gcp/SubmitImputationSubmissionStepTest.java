package bio.terra.pipelines.stairway.steps.imputation.gcp;

import static bio.terra.pipelines.testutils.TestUtils.VALID_METHOD_CONFIGURATION;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
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

class SubmitImputationSubmissionStepTest extends BaseEmbeddedDbTest {
  @Mock private RawlsService rawlsService;
  @Captor private ArgumentCaptor<SubmissionRequest> submissionRequestCaptor;
  @Captor private ArgumentCaptor<MethodConfiguration> updateMethodConfigCaptor;
  @Captor private ArgumentCaptor<MethodConfiguration> setMethodConfigCaptor;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;
  @Autowired private ImputationConfiguration imputationConfiguration;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final UUID randomUUID = UUID.randomUUID();

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    MethodConfiguration returnedMethodConfiguration = new MethodConfiguration();
    when(rawlsService.getCurrentMethodConfigForMethod(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.TEST_TOOL_NAME_1))
        .thenReturn(returnedMethodConfiguration);
    when(rawlsService.validateMethodConfig(
            returnedMethodConfiguration,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_TOOL_NAME_1,
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST,
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST,
            TestUtils.TEST_TOOL_VERSION_1))
        .thenReturn(true);
    when(rawlsService.submitWorkflow(
            eq("thisToken"),
            submissionRequestCaptor.capture(),
            eq(TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT),
            eq(TestUtils.CONTROL_WORKSPACE_NAME)))
        .thenReturn(new SubmissionReport().submissionId(randomUUID.toString()));

    // do the step
    SubmitImputationSubmissionStep submitImputationSubmissionStep =
        new SubmitImputationSubmissionStep(rawlsService, samService, imputationConfiguration);
    StepResult result = submitImputationSubmissionStep.doStep(flightContext);

    // extract the captured RunSetRequest and validate
    SubmissionRequest submissionRequest = submissionRequestCaptor.getValue();
    assertFalse(submissionRequest.isDeleteIntermediateOutputFiles());
    assertFalse(submissionRequest.isUseReferenceDisks());
    assertTrue(submissionRequest.isUseCallCache());
    assertEquals(testJobId.toString(), submissionRequest.getEntityName());

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        randomUUID,
        flightContext.getWorkingMap().get(ImputationJobMapKeys.SUBMISSION_ID, UUID.class));
  }

  @Test
  void doStepWithInvalidMethodConfig() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    // set up "current" method config to be version 1.1.1 with the corresponding method uri
    MethodConfiguration returnedMethodConfiguration =
        new MethodConfiguration()
            .methodRepoMethod(
                new MethodRepoMethod().methodUri("http/path/to/wdl/1.1.1").methodVersion("1.1.1"));
    when(rawlsService.getCurrentMethodConfigForMethod(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.TEST_TOOL_NAME_1))
        .thenReturn(returnedMethodConfiguration);
    when(rawlsService.validateMethodConfig(
            returnedMethodConfiguration,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_TOOL_NAME_1,
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST,
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST,
            TestUtils.TEST_TOOL_VERSION_1))
        .thenReturn(false);
    when(rawlsService.updateMethodConfigToBeValid(
            updateMethodConfigCaptor.capture(),
            eq(PipelinesEnum.ARRAY_IMPUTATION.getValue()),
            eq(TestUtils.TEST_TOOL_NAME_1),
            eq(TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST),
            eq(TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST),
            eq(TestUtils.TEST_TOOL_VERSION_1)))
        .thenReturn(VALID_METHOD_CONFIGURATION);
    when(rawlsService.setMethodConfigForMethod(
            eq("thisToken"),
            setMethodConfigCaptor.capture(),
            eq(TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT),
            eq(TestUtils.CONTROL_WORKSPACE_NAME),
            eq(TestUtils.TEST_TOOL_NAME_1)))
        .thenReturn(null);
    when(rawlsService.submitWorkflow(
            eq("thisToken"),
            any(SubmissionRequest.class),
            eq(TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT),
            eq(TestUtils.CONTROL_WORKSPACE_NAME)))
        .thenReturn(new SubmissionReport().submissionId(randomUUID.toString()));

    // do the step
    SubmitImputationSubmissionStep submitImputationSubmissionStep =
        new SubmitImputationSubmissionStep(rawlsService, samService, imputationConfiguration);
    submitImputationSubmissionStep.doStep(flightContext);

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
    MethodConfiguration returnedMethodConfiguration = new MethodConfiguration();
    when(rawlsService.getCurrentMethodConfigForMethod(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.TEST_TOOL_NAME_1))
        .thenReturn(returnedMethodConfiguration);
    when(rawlsService.validateMethodConfig(
            returnedMethodConfiguration,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_TOOL_NAME_1,
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST,
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST,
            TestUtils.TEST_TOOL_VERSION_1))
        .thenReturn(true);

    // throw exception on submitting workflow
    when(rawlsService.submitWorkflow(
            eq("thisToken"),
            any(SubmissionRequest.class),
            eq(TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT),
            eq(TestUtils.CONTROL_WORKSPACE_NAME)))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));
    // do the step
    SubmitImputationSubmissionStep submitImputationSubmissionStep =
        new SubmitImputationSubmissionStep(rawlsService, samService, imputationConfiguration);
    StepResult result = submitImputationSubmissionStep.doStep(flightContext);
    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());

    // throw exception on setting method config to be valid
    when(rawlsService.validateMethodConfig(
            returnedMethodConfiguration,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_TOOL_NAME_1,
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST,
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST,
            TestUtils.TEST_TOOL_VERSION_1))
        .thenReturn(false);
    when(rawlsService.setMethodConfigForMethod(
            "thisToken",
            null,
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.TEST_TOOL_NAME_1))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));
    // do the step
    submitImputationSubmissionStep =
        new SubmitImputationSubmissionStep(rawlsService, samService, imputationConfiguration);
    result = submitImputationSubmissionStep.doStep(flightContext);
    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());

    // throw exception getting method config
    when(rawlsService.getCurrentMethodConfigForMethod(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.TEST_TOOL_NAME_1))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));
    // do the step
    submitImputationSubmissionStep =
        new SubmitImputationSubmissionStep(rawlsService, samService, imputationConfiguration);
    result = submitImputationSubmissionStep.doStep(flightContext);
    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    SubmitImputationSubmissionStep submitImputationSubmissionStep =
        new SubmitImputationSubmissionStep(rawlsService, samService, imputationConfiguration);
    StepResult result = submitImputationSubmissionStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
