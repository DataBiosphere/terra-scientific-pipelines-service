package bio.terra.pipelines.stairway.steps.datadelivery;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

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

class UpdateDataDeliveryStatusStepTest extends BaseEmbeddedDbTest {

  @Mock private PipelineRunsService pipelineRunsService;
  @Mock private DataDeliveryService dataDeliveryService;
  @Mock private FlightContext flightContext;

  private final UUID testPipelineRunId = TestUtils.TEST_NEW_UUID;
  private final String testUserId = TestUtils.TEST_USER_1_ID;
  private PipelineRun testPipelineRun;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    inputParameters.put(JobMapKeys.USER_ID, testUserId);
    inputParameters.put(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, testPipelineRunId);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);

    testPipelineRun = TestUtils.createNewPipelineRunWithJobIdAndUser(testPipelineRunId, testUserId);
  }

  @Test
  void doStepSuccess() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId))
        .thenReturn(testPipelineRun);

    UpdateDataDeliveryStatusStep step =
        new UpdateDataDeliveryStatusStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
    verify(dataDeliveryService)
        .updateDataDeliveryStatus(testPipelineRun.getId(), DataDeliveryStatusEnum.SUCCEEDED);
  }

  @Test
  void doStepPipelineRunNotFound() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId)).thenReturn(null);

    UpdateDataDeliveryStatusStep step =
        new UpdateDataDeliveryStatusStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
  }

  @Test
  void undoStepSuccess() {
    UpdateDataDeliveryStatusStep step =
        new UpdateDataDeliveryStatusStep(pipelineRunsService, dataDeliveryService);
    StepResult result = step.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
