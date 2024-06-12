package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
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

    workingMap.put(RunImputationJobFlightMapKeys.WDS_URI, "wdsUri");

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() throws WdsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    Map<String, String> expectedOutputs =
        Map.of(
            "imputed_multi_sample_vcf", "some/file.vcf.gz",
            "imputed_multi_sample_vcf_index", "some/file.vcf.gz.tbi",
            "chunks_info", "some/file");
    RecordAttributes recordAttributes = new RecordAttributes();
    recordAttributes.putAll(expectedOutputs);
    RecordResponse recordResponse = new RecordResponse().attributes(recordAttributes);

    when(wdsService.getRecord(any(), any(), any(), any(), any())).thenReturn(recordResponse);

    FetchOutputsFromWdsStep fetchOutputsFromWdsStep =
        new FetchOutputsFromWdsStep(wdsService, samService);
    StepResult result = fetchOutputsFromWdsStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        expectedOutputs,
        flightContext
            .getWorkingMap()
            .get(RunImputationJobFlightMapKeys.RAW_OUTPUTS_MAP, Map.class));
  }

  @Test
  void doStepFailureRetry() throws WdsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    when(wdsService.getRecord(any(), any(), any(), any(), any()))
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
