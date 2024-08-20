package bio.terra.pipelines.stairway.imputation.gcp;

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
import bio.terra.rawls.model.SubmissionReport;
import bio.terra.rawls.model.SubmissionRequest;
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

class SubmitCromwellRunSetStepTest extends BaseEmbeddedDbTest {
  @Mock private RawlsService rawlsService;
  @Captor private ArgumentCaptor<SubmissionRequest> submissionRequestCaptor;
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
    when(rawlsService.submitWorkflow(any(), submissionRequestCaptor.capture(), any(), any()))
        .thenReturn(new SubmissionReport().submissionId(testJobId.toString()));

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(rawlsService, samService);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

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
  void doStepRawlsErrorRetry() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(rawlsService.submitWorkflow(any(), any(), any(), any()))
        .thenThrow(new RawlsServiceApiException("rawls is bad"));

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(rawlsService, samService);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // assert step is marked as retryable
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(rawlsService, samService);
    StepResult result = submitCromwellRunSetStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
