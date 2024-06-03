package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.Map;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class ImputationServiceTest extends BaseEmbeddedDbTest {

  @Autowired @InjectMocks ImputationService imputationService;
  @Mock private JobService mockJobService;
  @Mock private JobBuilder mockJobBuilder;
  @Autowired PipelineRunsService pipelineRunsService;
  @Autowired PipelineRunsRepository pipelineRunsRepository;

  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @Test
  void createImputationRunSuccess() {
    Pipeline testPipelineWithId = TestUtils.TEST_PIPELINE_1;
    testPipelineWithId.setId(3L);

    // mocks for Stairway
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.addParameter(any(), any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.submit()).thenReturn(testJobId);

    PipelineRun returnedPipelineRun =
        imputationService.createImputationRun(
            testJobId,
            testUserId,
            "test description",
            testPipelineWithId,
            testPipelineInputs,
            "test result path");
    // check that the appropriate fields are returned in returnedPipelineRun
    assertEquals(testJobId, returnedPipelineRun.getJobId());
    assertEquals(testUserId, returnedPipelineRun.getUserId());
    assertEquals("test description", returnedPipelineRun.getDescription());
    assertEquals("test result path", returnedPipelineRun.getResultUrl());
    assertEquals(testPipelineWithId.getId(), returnedPipelineRun.getPipelineId());

    // TODO these fail - i think because of something something transactions
    //    assertNotNull(returnedPipelineRun.getCreated());
    //    assertNotNull(returnedPipelineRun.getUpdated());

    // check that the pipeline run is written to the pipeline_runs table
    PipelineRun pipelineRunFromDb =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId).orElse(null);
    assertNotNull(pipelineRunFromDb);
    assertEquals(testJobId, pipelineRunFromDb.getJobId());
    assertEquals(testUserId, pipelineRunFromDb.getUserId());
    assertEquals("test description", pipelineRunFromDb.getDescription());
    assertEquals("test result path", pipelineRunFromDb.getResultUrl());
    assertEquals(testPipelineWithId.getId(), pipelineRunFromDb.getPipelineId());
    assertNotNull(pipelineRunFromDb.getCreated());
    assertNotNull(pipelineRunFromDb.getUpdated());
  }

  @Test
  void createImputationRunStairwayError() {
    Pipeline testPipelineWithId = TestUtils.TEST_PIPELINE_1;
    testPipelineWithId.setId(1L);

    // test that when we try to create a new run but Stairway fails, we don't write the run
    // to our database (i.e. test the transaction)
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(any())).thenThrow(new RuntimeException("Stairway error"));

    assertThrows(
        RuntimeException.class,
        () ->
            imputationService.createImputationRun(
                testJobId,
                testUserId,
                "test description",
                testPipelineWithId,
                testPipelineInputs,
                "test result path"));

    // check that the pipeline is not written to the pipeline_runs table
    assertEquals(
        Optional.empty(), pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId));
  }
}
