package bio.terra.pipelines.stairway.imputation.steps.azure;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
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
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() {
    // setup
    List<ListAppResponse> returnedListAppListResponse =
        List.of(new ListAppResponse().appName("justSomething"));
    when(leonardoService.getApps(TestUtils.CONTROL_WORKSPACE_ID.toString(), "thisToken"))
        .thenReturn(returnedListAppListResponse);
    when(leonardoService.getCbasUrlFromGetAppResponse(
            returnedListAppListResponse, TestUtils.CONTROL_WORKSPACE_ID.toString()))
        .thenReturn("cbasUriRetrieved");
    when(leonardoService.getWdsUrlFromGetAppResponse(
            returnedListAppListResponse, TestUtils.CONTROL_WORKSPACE_ID.toString()))
        .thenReturn("wdsUriRetrieved");

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
  void undoStepSuccess() {
    GetAppUrisStep getAppUrisStep = new GetAppUrisStep(leonardoService, samService);
    StepResult result = getAppUrisStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
