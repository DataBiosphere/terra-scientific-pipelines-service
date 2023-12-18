package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class JobsServiceMockTest extends BaseContainerTest {

  @Mock private StairwayJobService mockStairwayJobService;
  @Mock private StairwayJobBuilder mockStairwayJobBuilder;

  @InjectMocks private JobsService jobsService;

  // parameters used repeatedly by various tests, and things we'll want mocks to respond to
  // universally
  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final String testGoodPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;

  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testUUID = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockStairwayJobService.newJob()).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.jobId(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.flightClass(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.pipelineId(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.pipelineVersion(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.submittingUserId(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.pipelineInputs(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.submit()).thenReturn(testUUID);
  }

  @Test
  void testCreateJob_success() {
    // a job isn't actually kicked off
    UUID writtenUUID =
        jobsService.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);
    assertEquals(testUUID, writtenUUID);
  }
}
