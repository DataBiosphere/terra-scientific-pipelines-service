package bio.terra.pipelines.stairway.steps.datadelivery;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.verifyNoInteractions;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.DataDeliveryStatusEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.DataDeliveryService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.flights.datadelivery.DataDeliveryJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class CreateDataDeliveryRecordStepTest extends BaseEmbeddedDbTest {

  @Mock private PipelineRunsService pipelineRunsService;
  @Mock private DataDeliveryService dataDeliveryService;
  @Mock private FlightContext flightContext;

  private final UUID testPipelineRunId = TestUtils.TEST_NEW_UUID;
  private final String testUserId = TestUtils.TEST_USER_1_ID;
  private final GcsFile testDestinationGcsPath = new GcsFile("gs://test-bucket/test-path");
  private PipelineRun testPipelineRun;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    inputParameters.put(JobMapKeys.USER_ID, testUserId);
    inputParameters.put(JobMapKeys.DOMAIN_NAME, TestUtils.TEST_DOMAIN);
    inputParameters.put(DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, testDestinationGcsPath);
    inputParameters.put(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, testPipelineRunId);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getFlightId()).thenReturn(testPipelineRunId.toString());

    testPipelineRun = TestUtils.createNewPipelineRunWithJobIdAndUser(testPipelineRunId, testUserId);
  }

  @Test
  void doStepSuccess() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId))
        .thenReturn(testPipelineRun);

    CreateDataDeliveryRecordStep step =
        new CreateDataDeliveryRecordStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
    verify(dataDeliveryService)
        .createDataDelivery(
            testPipelineRun.getId(),
            testPipelineRunId,
            DataDeliveryStatusEnum.RUNNING,
            testDestinationGcsPath);
  }

  @Test
  void doStepPipelineRunNotFound() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId)).thenReturn(null);

    CreateDataDeliveryRecordStep step =
        new CreateDataDeliveryRecordStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
    verifyNoInteractions(dataDeliveryService);
  }

  @Test
  void undoStepPipelineRunFound() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId))
        .thenReturn(testPipelineRun);

    CreateDataDeliveryRecordStep step =
        new CreateDataDeliveryRecordStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
    verify(dataDeliveryService)
        .updateDataDeliveryStatus(testPipelineRun.getId(), DataDeliveryStatusEnum.FAILED);
  }

  @Test
  void undoStepPipelineRunNotFound() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId)).thenReturn(null);

    CreateDataDeliveryRecordStep step =
        new CreateDataDeliveryRecordStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verifyNoInteractions(dataDeliveryService);
  }
}
