package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class ImputationServiceMockTest extends BaseEmbeddedDbTest {
  @InjectMocks private ImputationService imputationService;
  @Mock private JobService mockJobService;
  @Mock private JobBuilder mockJobBuilder;

  // parameters used repeatedly by various tests, and things we'll want mocks to respond to
  // universally
  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testUUID = TestUtils.TEST_NEW_UUID;

  @Mock ImputationConfiguration imputationConfiguration;

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.addParameter(any(), any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.submit()).thenReturn(testUUID);

    imputationConfiguration.setCromwellSubmissionPollingIntervalInSeconds(1L);
  }

  @Test
  void createJobSuccess() {
    // note this doesn't actually kick off a job
    UUID writtenUUID =
        imputationService.createImputationJob(
            testUUID,
            testUserId,
            "test description",
            TestUtils.TEST_PIPELINE_1,
            testPipelineInputs,
            TestUtils.TEST_RESULT_URL);
    assertEquals(testUUID, writtenUUID);
  }
}
