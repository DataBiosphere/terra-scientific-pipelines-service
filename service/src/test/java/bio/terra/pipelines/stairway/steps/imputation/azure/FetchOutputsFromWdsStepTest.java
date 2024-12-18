package bio.terra.pipelines.stairway.steps.imputation.azure;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.Map;
import org.databiosphere.workspacedata.client.ApiException;
import org.databiosphere.workspacedata.model.RecordAttributes;
import org.databiosphere.workspacedata.model.RecordResponse;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class FetchOutputsFromWdsStepTest extends BaseEmbeddedDbTest {

  @Mock WdsService wdsService;
  @Mock SamService samService;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    workingMap.put(ImputationJobMapKeys.WDS_URI, "wdsUri");

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() throws WdsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // outputs to match the test output definitions
    Map<String, String> expectedOutputsFromWds = Map.of("output_name", "some/file.vcf.gz");
    RecordAttributes recordAttributes = new RecordAttributes();
    recordAttributes.putAll(expectedOutputsFromWds);
    RecordResponse recordResponse = new RecordResponse().attributes(recordAttributes);

    when(wdsService.getRecord(
            "wdsUri",
            "thisToken",
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.CONTROL_WORKSPACE_ID.toString(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenReturn(recordResponse);

    FetchOutputsFromWdsStep fetchOutputsFromWdsStep =
        new FetchOutputsFromWdsStep(wdsService, samService);
    StepResult result = fetchOutputsFromWdsStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // in the step we translate the snake_case output keys to camelCase
    Map<String, String> expectedOutputsFromWorkingMap = Map.of("outputName", "some/file.vcf.gz");

    assertEquals(
        expectedOutputsFromWorkingMap,
        flightContext.getWorkingMap().get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class));
  }

  @Test
  void doStepFailureRetry() throws WdsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    when(wdsService.getRecord(
            "wdsUri",
            "thisToken",
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.CONTROL_WORKSPACE_ID.toString(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenThrow(new WdsServiceApiException(new ApiException("WDS Service Api Exception")));

    FetchOutputsFromWdsStep fetchOutputsFromWdsStep =
        new FetchOutputsFromWdsStep(wdsService, samService);
    StepResult result = fetchOutputsFromWdsStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    FetchOutputsFromWdsStep fetchOutputsFromWdsStep =
        new FetchOutputsFromWdsStep(wdsService, samService);
    StepResult result = fetchOutputsFromWdsStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
