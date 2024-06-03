package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.time.LocalDateTime;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class ImputationServiceMockTest extends BaseEmbeddedDbTest {
  @InjectMocks private ImputationService imputationService;
  @Mock private PipelineRunsService mockPipelineRunsService;
  @Mock private JobService mockJobService;
  @Mock private JobBuilder mockJobBuilder;

  // parameters used repeatedly by various tests, and things we'll want mocks to respond to
  // universally
  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testUUID = TestUtils.TEST_NEW_UUID;
  private final LocalDateTime createdTime = LocalDateTime.now();
  private final LocalDateTime updatedTime = LocalDateTime.now();
  private final PipelineRun expectedPipelineRun =
      new PipelineRun(
          testUUID,
          testUserId,
          1L,
          createdTime,
          updatedTime,
          CommonPipelineRunStatusEnum.RUNNING.toString(),
          "test description",
          TestUtils.TEST_RESULT_URL,
          null,
          null);

  @Mock ImputationConfiguration imputationConfiguration;

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.addParameter(any(), any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.submit()).thenReturn(testUUID);

    when(mockPipelineRunsService.writePipelineRunToDb(
            eq(testUUID),
            any(String.class),
            any(Long.class),
            any(CommonPipelineRunStatusEnum.class),
            any(String.class),
            any(String.class),
            any(Map.class)))
        .thenReturn(expectedPipelineRun);

    imputationConfiguration.setCromwellSubmissionPollingIntervalInSeconds(1L);
  }

  @Test
  void createImputationRunSuccess() {
    Pipeline testPipelineWithId = TestUtils.TEST_PIPELINE_1;
    testPipelineWithId.setId(1L);
    // note this test doesn't kick off a real job
    PipelineRun writtenPipelineRun =
        imputationService.createImputationRun(
            testUUID,
            testUserId,
            "test description",
            testPipelineWithId,
            testPipelineInputs,
            TestUtils.TEST_RESULT_URL);

    assertEquals(expectedPipelineRun, writtenPipelineRun);
  }

  @Test
  void createImputationRunDbWriteFail() {
    Pipeline testPipelineWithId = TestUtils.TEST_PIPELINE_1;
    testPipelineWithId.setId(1L);
    // override this mock to throw an error
    when(mockPipelineRunsService.writePipelineRunToDb(
            eq(testUUID),
            any(String.class),
            any(Long.class),
            any(CommonPipelineRunStatusEnum.class),
            any(String.class),
            any(String.class),
            any(Map.class)))
        .thenThrow(new RuntimeException("test exception"));

    assertThrows(
        RuntimeException.class,
        () ->
            imputationService.createImputationRun(
                testUUID,
                testUserId,
                "test description",
                testPipelineWithId,
                testPipelineInputs,
                TestUtils.TEST_RESULT_URL));
  }

  @Test
  void createImputationRunStairwayError() throws Exception {
    Pipeline testPipelineWithId = TestUtils.TEST_PIPELINE_1;
    testPipelineWithId.setId(1L);
    // override this mock to throw an error
    when(mockJobBuilder.submit()).thenThrow(new RuntimeException("test exception"));

    assertThrows(
        RuntimeException.class,
        () ->
            imputationService.createImputationRun(
                testUUID,
                testUserId,
                "test description",
                testPipelineWithId,
                testPipelineInputs,
                TestUtils.TEST_RESULT_URL));
  }
}
