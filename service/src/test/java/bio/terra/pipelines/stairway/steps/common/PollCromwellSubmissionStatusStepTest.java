package bio.terra.pipelines.stairway.steps.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.Submission;
import bio.terra.rawls.model.SubmissionStatus;
import bio.terra.rawls.model.Workflow;
import bio.terra.rawls.model.WorkflowStatus;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class PollCromwellSubmissionStatusStepTest extends BaseEmbeddedDbTest {

  @Mock private RawlsService rawlsService;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final UUID randomUUID = UUID.randomUUID();
  private final String toolConfigKey = TestUtils.TOOL_CONFIG_KEY;
  private final String submissionIdKey = TestUtils.SUBMISSION_ID_KEY;
  private final ToolConfig toolConfig = TestUtils.TOOL_CONFIG_GENERIC;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    inputParameters.put(toolConfigKey, toolConfig);
    FlightMap workingMap = new FlightMap();
    workingMap.put(submissionIdKey, randomUUID);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    Submission response =
        new Submission()
            .status(SubmissionStatus.DONE)
            .addWorkflowsItem(new Workflow().status(WorkflowStatus.SUCCEEDED));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.getSubmissionStatus(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            randomUUID))
        .thenReturn(response);

    // do the step
    PollCromwellSubmissionStatusStep pollCromwellSubmissionStatusStep =
        new PollCromwellSubmissionStatusStep(
            rawlsService, samService, toolConfigKey, submissionIdKey);
    StepResult result = pollCromwellSubmissionStatusStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepRunningThenComplete() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    Submission firstResponse =
        new Submission()
            .status(SubmissionStatus.SUBMITTED)
            .addWorkflowsItem(new Workflow().status(WorkflowStatus.RUNNING));
    Submission secondResponse =
        new Submission()
            .status(SubmissionStatus.DONE)
            .addWorkflowsItem(new Workflow().status(WorkflowStatus.SUCCEEDED));

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.getSubmissionStatus(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            randomUUID))
        .thenReturn(firstResponse)
        .thenReturn(secondResponse);

    // do the step
    PollCromwellSubmissionStatusStep pollCromwellSubmissionStatusStep =
        new PollCromwellSubmissionStatusStep(
            rawlsService, samService, toolConfigKey, submissionIdKey);
    StepResult result = pollCromwellSubmissionStatusStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepNotAllSuccessfulRuns() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    Submission responseWithErrorRun =
        new Submission()
            .status(SubmissionStatus.DONE)
            .addWorkflowsItem(new Workflow().status(WorkflowStatus.SUCCEEDED))
            .addWorkflowsItem(new Workflow().status(WorkflowStatus.FAILED));

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.getSubmissionStatus(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            randomUUID))
        .thenReturn(responseWithErrorRun);

    // do the step
    PollCromwellSubmissionStatusStep pollCromwellSubmissionStatusStep =
        new PollCromwellSubmissionStatusStep(
            rawlsService, samService, toolConfigKey, submissionIdKey);
    StepResult result = pollCromwellSubmissionStatusStep.doStep(flightContext);

    // make sure the step fails and the exception message contains relevant info
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "An error occurred while running the job. Not all workflows succeeded in submission. Please contact support for help.",
        result.getException().get().getMessage());
  }

  @Test
  void doStepRawlsApiErrorRetry() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.getSubmissionStatus(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            randomUUID))
        .thenThrow(new RawlsServiceApiException("this is the error message"));

    // do the step, expect a Retry status
    PollCromwellSubmissionStatusStep pollCromwellSubmissionStatusStep =
        new PollCromwellSubmissionStatusStep(
            rawlsService, samService, toolConfigKey, submissionIdKey);
    StepResult result = pollCromwellSubmissionStatusStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    PollCromwellSubmissionStatusStep pollCromwellSubmissionStatusStep =
        new PollCromwellSubmissionStatusStep(
            rawlsService, samService, toolConfigKey, submissionIdKey);
    StepResult result = pollCromwellSubmissionStatusStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
