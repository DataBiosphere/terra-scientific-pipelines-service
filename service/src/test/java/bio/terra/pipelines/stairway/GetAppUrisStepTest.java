package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.GetAppUrisStep;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
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

class GetAppUrisStepTest extends BaseEmbeddedDbTest {

  @Mock private LeonardoService leonardoService;
  @Mock private SamService samService;
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
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(leonardoService.checkHealth())
        .thenReturn(new HealthCheck.Result(true, "leonardo is healthy"));
    when(leonardoService.getCbasUrlFromGetAppResponse(any(), any())).thenReturn("cbasUriRetrieved");
    when(leonardoService.getWdsUrlFromGetAppResponse(any(), any())).thenReturn("wdsUriRetrieved");

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    GetAppUrisStep getAppUrisStep = new GetAppUrisStep(leonardoService, samService);
    StepResult result = getAppUrisStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap workingMap = flightContext.getWorkingMap();

    // make sure the working map was updated appropriately
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        "wdsUriRetrieved", workingMap.get(RunImputationJobFlightMapKeys.WDS_URI, String.class));
    assertEquals(
        "cbasUriRetrieved", workingMap.get(RunImputationJobFlightMapKeys.CBAS_URI, String.class));
  }

  @Test
  void doStepUnhealthyLeonardo() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(leonardoService.checkHealth())
        .thenReturn(new HealthCheck.Result(false, "leonardo is not healthy"));

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    GetAppUrisStep getAppUrisStep = new GetAppUrisStep(leonardoService, samService);
    StepResult result = getAppUrisStep.doStep(flightContext);

    // make sure the appropriate step status was returned
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() throws InterruptedException {
    GetAppUrisStep getAppUrisStep = new GetAppUrisStep(leonardoService, samService);
    StepResult result = getAppUrisStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
