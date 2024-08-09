package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class CheckWdsHealthStepTest extends BaseEmbeddedDbTest {
  @Mock private WdsService wdsService;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();
    workingMap.put(RunImputationJobFlightMapKeys.WDS_URI, "wdsUri");

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() {
    // setup
    when(wdsService.checkHealth(any(), any()))
        .thenReturn(new HealthCheckWorkspaceApps.Result(true, "wds is healthy"));

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    CheckWdsHealthStep checkWdsHealthStep = new CheckWdsHealthStep(wdsService, samService);
    StepResult result = checkWdsHealthStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepUnhealthyWds() {
    // setup
    when(wdsService.checkHealth(any(), any()))
        .thenReturn(new HealthCheckWorkspaceApps.Result(false, "wds is not healthy"));

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    CheckWdsHealthStep checkWdsHealthStep = new CheckWdsHealthStep(wdsService, samService);
    StepResult result = checkWdsHealthStep.doStep(flightContext);

    // make sure the appropriate step status was returned
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    CheckWdsHealthStep checkWdsHealthStep = new CheckWdsHealthStep(wdsService, samService);
    StepResult result = checkWdsHealthStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
