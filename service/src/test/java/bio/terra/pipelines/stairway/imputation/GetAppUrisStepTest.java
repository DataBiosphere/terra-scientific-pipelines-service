package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class GetAppUrisStepTest extends BaseEmbeddedDbTest {

  @Mock private LeonardoService leonardoService;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;

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
        "wdsUriRetrieved",
        workingMap.get(RunImputationAzureJobFlightMapKeys.WDS_URI, String.class));
    assertEquals(
        "cbasUriRetrieved",
        workingMap.get(RunImputationAzureJobFlightMapKeys.CBAS_URI, String.class));
  }

  @Test
  void undoStepSuccess() {
    GetAppUrisStep getAppUrisStep = new GetAppUrisStep(leonardoService, samService);
    StepResult result = getAppUrisStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
