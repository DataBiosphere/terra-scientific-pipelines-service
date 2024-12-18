package bio.terra.pipelines.stairway.steps.imputation.azure;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
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

class CheckLeonardoHealthStepTest extends BaseEmbeddedDbTest {

  @Mock private LeonardoService leonardoService;
  @Mock private FlightContext flightContext;

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
    when(leonardoService.checkHealth())
        .thenReturn(new HealthCheck.Result(true, "leonardo is healthy"));

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    CheckLeonardoHealthStep checkLeonardoHealthStep = new CheckLeonardoHealthStep(leonardoService);
    StepResult result = checkLeonardoHealthStep.doStep(flightContext);

    // make sure the working map was updated appropriately
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepUnhealthyLeonardo() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(leonardoService.checkHealth())
        .thenReturn(new HealthCheck.Result(false, "leonardo is not healthy"));

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    CheckLeonardoHealthStep checkLeonardoHealthStep = new CheckLeonardoHealthStep(leonardoService);
    StepResult result = checkLeonardoHealthStep.doStep(flightContext);

    // make sure the appropriate step status was returned
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    CheckLeonardoHealthStep checkLeonardoHealthStep = new CheckLeonardoHealthStep(leonardoService);
    StepResult result = checkLeonardoHealthStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
