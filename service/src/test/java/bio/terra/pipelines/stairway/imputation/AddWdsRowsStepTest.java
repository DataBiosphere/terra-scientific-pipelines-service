package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.*;
import java.util.UUID;
import org.databiosphere.workspacedata.client.ApiException;
import org.databiosphere.workspacedata.model.RecordAttributes;
import org.databiosphere.workspacedata.model.RecordRequest;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Captor;
import org.mockito.Mock;

class AddWdsRowsStepTest extends BaseEmbeddedDbTest {
  @Mock private WdsService wdsService;
  @Captor ArgumentCaptor<RecordRequest> recordRequestCaptor;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    workingMap.put(RunImputationJobFlightMapKeys.WDS_URI, "wdsUri");
    workingMap.put(
        RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, TestUtils.TEST_PIPELINE_INPUTS);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() throws WdsServiceException, InterruptedException {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    AddWdsRowStep addWdsRowStep = new AddWdsRowStep(wdsService, samService);
    StepResult result = addWdsRowStep.doStep(flightContext);

    // extract the captured RecordRequest
    verify(wdsService)
        .createOrReplaceRecord(
            any(), any(), recordRequestCaptor.capture(), any(), any(), any(), any());
    RecordRequest recordRequest = recordRequestCaptor.getValue();
    RecordAttributes recordAttributes = recordRequest.getAttributes();
    // record request should have all pipeline inputs plus a timestamp field
    assertEquals(TestUtils.TEST_PIPELINE_INPUTS.size() + 1, recordAttributes.size());
    for (String key : TestUtils.TEST_PIPELINE_INPUTS.keySet()) {
      assertEquals(TestUtils.TEST_PIPELINE_INPUTS.get(key), recordAttributes.get(key));
    }

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepWdsException() throws WdsServiceException, InterruptedException {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    WdsServiceApiException thrownException =
        new WdsServiceApiException(new ApiException("this is the error message"));
    when(wdsService.createOrReplaceRecord(any(), any(), any(), any(), any(), any(), any()))
        .thenThrow(thrownException);

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    AddWdsRowStep addWdsRowStep = new AddWdsRowStep(wdsService, samService);
    StepResult result = addWdsRowStep.doStep(flightContext);

    // make sure the step was fails with an exception
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(thrownException, result.getException().get());
  }

  @Test
  void undoStepSuccess() throws InterruptedException {
    AddWdsRowStep addWdsRowStep = new AddWdsRowStep(wdsService, samService);
    StepResult result = addWdsRowStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
