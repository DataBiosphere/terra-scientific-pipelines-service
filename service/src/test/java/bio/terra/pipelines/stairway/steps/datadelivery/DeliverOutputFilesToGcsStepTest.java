package bio.terra.pipelines.stairway.steps.datadelivery;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
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

class DeliverOutputFilesToGcsStepTest extends BaseEmbeddedDbTest {

  @Mock private PipelineRunsService pipelineRunsService;
  @Mock private PipelineInputsOutputsService pipelineInputsOutputsService;
  @Mock private FlightContext flightContext;

  private final UUID testPipelineRunId = TestUtils.TEST_NEW_UUID;
  private final String testUserId = TestUtils.TEST_USER_1_ID;
  private final GcsFile testDestinationGcsPath = new GcsFile("gs://test-bucket/test-path");
  private PipelineRun testPipelineRun;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    inputParameters.put(JobMapKeys.USER_ID, testUserId);
    inputParameters.put(JobMapKeys.PIPELINE_ID, TestUtils.TEST_PIPELINE_ID_1);
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

    DeliverOutputFilesToGcsStep step =
        new DeliverOutputFilesToGcsStep(pipelineRunsService, pipelineInputsOutputsService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
    verify(pipelineInputsOutputsService)
        .deliverOutputFilesToGcs(testPipelineRun, testDestinationGcsPath);
  }

  @Test
  void doStepPipelineRunNotFound() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId)).thenReturn(null);

    DeliverOutputFilesToGcsStep step =
        new DeliverOutputFilesToGcsStep(pipelineRunsService, pipelineInputsOutputsService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
  }

  @Test
  void doStepDeliveryFailureRetry() {
    when(pipelineRunsService.getPipelineRun(testPipelineRunId, testUserId))
        .thenReturn(testPipelineRun);
    doThrow(new RuntimeException("could not deliver files to destination"))
        .when(pipelineInputsOutputsService)
        .deliverOutputFilesToGcs(testPipelineRun, testDestinationGcsPath);

    DeliverOutputFilesToGcsStep step =
        new DeliverOutputFilesToGcsStep(pipelineRunsService, pipelineInputsOutputsService);
    StepResult result = step.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
    verify(pipelineRunsService).getPipelineRun(testPipelineRunId, testUserId);
    verify(pipelineInputsOutputsService)
        .deliverOutputFilesToGcs(testPipelineRun, testDestinationGcsPath);
  }

  @Test
  void undoStepSuccess() {
    DeliverOutputFilesToGcsStep step =
        new DeliverOutputFilesToGcsStep(pipelineRunsService, pipelineInputsOutputsService);
    StepResult result = step.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
