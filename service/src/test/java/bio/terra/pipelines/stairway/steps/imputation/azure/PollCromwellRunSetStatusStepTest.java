package bio.terra.pipelines.stairway.steps.imputation.azure;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.*;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class PollCromwellRunSetStatusStepTest extends BaseEmbeddedDbTest {

  @Mock private CbasService cbasService;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;
  @Autowired ImputationConfiguration imputationConfiguration;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final UUID randomUUID = UUID.randomUUID();

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();
    workingMap.put(ImputationJobMapKeys.CBAS_URI, "cbasUri");
    workingMap.put(ImputationJobMapKeys.RUN_SET_ID, randomUUID);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    RunLogResponse response =
        new RunLogResponse().addRunsItem(new RunLog().state(RunState.COMPLETE));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.getRunsForRunSet("cbasUri", "thisToken", randomUUID)).thenReturn(response);

    // do the step
    PollCromwellRunSetStatusStep pollCromwellRunSetStatusStep =
        new PollCromwellRunSetStatusStep(cbasService, samService, imputationConfiguration);
    StepResult result = pollCromwellRunSetStatusStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepRunningThenComplete() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    RunLogResponse firstResponse =
        new RunLogResponse()
            .addRunsItem(new RunLog().state(RunState.RUNNING))
            .addRunsItem(new RunLog().state(RunState.COMPLETE));
    RunLogResponse secondResponse =
        new RunLogResponse()
            .addRunsItem(new RunLog().state(RunState.COMPLETE))
            .addRunsItem(new RunLog().state(RunState.COMPLETE));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.checkHealth("cbasUri", "thisToken"))
        .thenReturn(new HealthCheckWorkspaceApps.Result(true, "cbas is healthy"));
    when(cbasService.getRunsForRunSet("cbasUri", "thisToken", randomUUID))
        .thenReturn(firstResponse)
        .thenReturn(secondResponse);

    // do the step
    PollCromwellRunSetStatusStep pollCromwellRunSetStatusStep =
        new PollCromwellRunSetStatusStep(cbasService, samService, imputationConfiguration);
    StepResult result = pollCromwellRunSetStatusStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepCbasApiErrorRetry() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.getRunsForRunSet("cbasUri", "thisToken", randomUUID))
        .thenThrow(new CbasServiceApiException("this is the error message"));

    // do the step, expect a Retry status
    PollCromwellRunSetStatusStep pollCromwellRunSetStatusStep =
        new PollCromwellRunSetStatusStep(cbasService, samService, imputationConfiguration);
    StepResult result = pollCromwellRunSetStatusStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void doStepNotAllSuccessfulRuns() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    RunLogResponse responseWithErrorRun =
        new RunLogResponse()
            .addRunsItem(new RunLog().state(RunState.COMPLETE))
            .addRunsItem(new RunLog().state(RunState.EXECUTOR_ERROR));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.getRunsForRunSet("cbasUri", "thisToken", randomUUID))
        .thenReturn(responseWithErrorRun);

    // do the step
    PollCromwellRunSetStatusStep pollCromwellRunSetStatusStep =
        new PollCromwellRunSetStatusStep(cbasService, samService, imputationConfiguration);
    StepResult result = pollCromwellRunSetStatusStep.doStep(flightContext);

    // make sure the step fails
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    PollCromwellRunSetStatusStep pollCromwellRunSetStatusStep =
        new PollCromwellRunSetStatusStep(cbasService, samService, imputationConfiguration);
    StepResult result = pollCromwellRunSetStatusStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
