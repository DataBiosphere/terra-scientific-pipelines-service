package bio.terra.pipelines.stairway.imputation.gcp;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.client.ApiException;
import bio.terra.rawls.model.Entity;
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

class AddDataTableRowStepTest extends BaseEmbeddedDbTest {
  @Mock private RawlsService rawlsService;
  @Captor ArgumentCaptor<Entity> entityCaptor;
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
  void doStepSuccess() throws RawlsServiceException {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    AddDataTableRowStep addDataTableRowStep = new AddDataTableRowStep(rawlsService, samService);
    StepResult result = addDataTableRowStep.doStep(flightContext);

    // extract the captured RecordRequest
    verify(rawlsService).upsertDataTableEntity(any(), any(), any(), entityCaptor.capture());
    Entity capturedEntity = entityCaptor.getValue();
    // validate fields of captured entity
    assertEquals(flightContext.getFlightId(), capturedEntity.getName());

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepRawlsExceptionRetry() throws RawlsServiceApiException {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    RawlsServiceApiException thrownException =
        new RawlsServiceApiException(new ApiException("this is the error message"));
    when(rawlsService.upsertDataTableEntity(any(), any(), any(), any())).thenThrow(thrownException);

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step, expect a Retry status
    AddDataTableRowStep addDataTableRowStep = new AddDataTableRowStep(rawlsService, samService);
    StepResult result = addDataTableRowStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    AddDataTableRowStep addDataTableRowStep = new AddDataTableRowStep(rawlsService, samService);
    StepResult result = addDataTableRowStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
