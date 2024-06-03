package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlight;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelineRunsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelineRunsService pipelineRunsService;
  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @Autowired PipelineInputsRepository pipelineInputsRepository;

  // mock Stairway
  @MockBean private JobService mockJobService;
  @MockBean private JobBuilder mockJobBuilder;

  private final String testUserId = TestUtils.TEST_USER_ID_1;

  private final CommonPipelineRunStatusEnum testStatus = CommonPipelineRunStatusEnum.RUNNING;
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final String testDescription = TestUtils.TEST_PIPELINE_DESCRIPTION_1;
  private final String testResultUrl = TestUtils.TEST_RESULT_URL;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private PipelineRun createTestJobWithJobId(UUID jobId) {
    return createTestJobWithJobIdAndUser(jobId, testUserId);
  }

  private PipelineRun createTestJobWithJobIdAndUser(UUID jobId, String userId) {
    return new PipelineRun(
        jobId, userId, testPipelineId, testStatus.toString(), testDescription, testResultUrl);
  }

  private Pipeline createTestPipelineWithId() {
    Pipeline pipeline =
        new Pipeline(
            TestUtils.TEST_PIPELINE_1.getName(),
            TestUtils.TEST_PIPELINE_1.getVersion(),
            TestUtils.TEST_PIPELINE_1.getDisplayName(),
            TestUtils.TEST_PIPELINE_1.getDescription(),
            TestUtils.TEST_PIPELINE_1.getPipelineType(),
            TestUtils.TEST_PIPELINE_1.getWdlUrl(),
            TestUtils.TEST_PIPELINE_1.getWdlMethodName(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceId(),
            TestUtils.TEST_PIPELINE_1.getPipelineInputDefinitions());
    pipeline.setId(3L);
    return pipeline;
  }

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.addParameter(any(), any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.submit()).thenReturn(testJobId);
  }

  @Test
  void writeJobToDbOk() {
    List<PipelineRun> runsDefault = pipelineRunsRepository.findAllByUserId(testUserId);

    // test data migration inserts one row by default
    assertEquals(1, runsDefault.size());

    PipelineRun pipelineRun =
        pipelineRunsService.writePipelineRunToDb(
            testJobId,
            testUserId,
            testPipelineId,
            testStatus,
            testDescription,
            testResultUrl,
            testPipelineInputs);

    List<PipelineRun> runsAfterSave = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(2, runsAfterSave.size());

    // verify info written to the pipeline_runs table
    PipelineRun savedRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(pipelineRun.getJobId(), testUserId)
            .orElseThrow();
    assertEquals(testJobId, savedRun.getJobId());
    assertEquals(testPipelineId, savedRun.getPipelineId());
    assertEquals(testUserId, savedRun.getUserId());
    assertEquals(testDescription, savedRun.getDescription());
    assertEquals(testResultUrl, savedRun.getResultUrl());
    assertEquals(testStatus.toString(), savedRun.getStatus());
    assertNotNull(savedRun.getCreated());
    assertNotNull(savedRun.getUpdated());

    // verify info written to pipelineInputs table
    Optional<PipelineInput> pipelineInput = pipelineInputsRepository.findById(savedRun.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals("{\"first_key\":\"first_value\"}", pipelineInput.get().getInputs());
  }

  @Test
  void writeRunToDbDuplicateRun() {
    // try to save a run with the same job id two times, the second time it should throw duplicate
    // exception error
    PipelineRun newPipelineRun = createTestJobWithJobId(testJobId);

    PipelineRun savedJobFirst =
        pipelineRunsService.writePipelineRunToDbThrowsDuplicateException(newPipelineRun);
    assertNotNull(savedJobFirst);

    PipelineRun newPipelineRunSameId = createTestJobWithJobId(testJobId);
    assertThrows(
        DuplicateObjectException.class,
        () ->
            pipelineRunsService.writePipelineRunToDbThrowsDuplicateException(newPipelineRunSameId));
  }

  @Test
  void testGetCorrectNumberOfRows() {
    // A test row should exist for this user.
    List<PipelineRun> pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(1, pipelineRuns.size());

    // insert another row and verify that it shows up
    PipelineRun newPipelineRun = createTestJobWithJobId(testJobId);

    pipelineRunsRepository.save(newPipelineRun);
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(2, pipelineRuns.size());
  }

  @Test
  void testCorrectUserIsolation() {
    // A test row should exist for this user.
    List<PipelineRun> pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(1, pipelineRuns.size());

    // insert row for second user and verify that it shows up
    String testUserId2 = TestUtils.TEST_USER_ID_2;
    PipelineRun newPipelineRun = createTestJobWithJobIdAndUser(UUID.randomUUID(), testUserId2);
    pipelineRunsRepository.save(newPipelineRun);

    // Verify that the old userid still show only 1 record
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(1, pipelineRuns.size());

    // Verify the new user's id shows a single job as well
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId2);
    assertEquals(1, pipelineRuns.size());
  }

  @Test
  void createPipelineRunNoWorkspaceSetUp() {
    Pipeline testPipelineWithIdMissingWorkspaceId = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceId.setWorkspaceId(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.createPipelineRun(
                testPipelineWithIdMissingWorkspaceId,
                testJobId,
                testUserId,
                testDescription,
                testPipelineInputs,
                testResultUrl));
  }

  @Test
  void createPipelineRunImputation() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // override this mock to ensure the correct flight class is being requested
    when(mockJobBuilder.flightClass(RunImputationJobFlight.class)).thenReturn(mockJobBuilder);

    PipelineRun returnedPipelineRun =
        pipelineRunsService.createPipelineRun(
            testPipelineWithId,
            testJobId,
            testUserId,
            testDescription,
            testPipelineInputs,
            testResultUrl);

    assertEquals(testJobId, returnedPipelineRun.getJobId());
    assertEquals(testUserId, returnedPipelineRun.getUserId());
    assertEquals(testDescription, returnedPipelineRun.getDescription());
    assertEquals(testResultUrl, returnedPipelineRun.getResultUrl());
    assertEquals(testPipelineWithId.getId(), returnedPipelineRun.getPipelineId());
    //    assertNotNull(returnedPipelineRun.getCreated());
    //    assertNotNull(returnedPipelineRun.getUpdated());

    // verify info written to pipeline_runs table
    PipelineRun savedRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(returnedPipelineRun.getJobId(), testUserId)
            .orElseThrow();
    assertEquals(testJobId, savedRun.getJobId());
    assertNotNull(savedRun.getCreated());
    assertNotNull(savedRun.getUpdated());

    // verify info written to pipeline_inputs table
    Optional<PipelineInput> pipelineInput =
        pipelineInputsRepository.findById(returnedPipelineRun.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals("{\"first_key\":\"first_value\"}", pipelineInput.get().getInputs());

    // verify submit was called
    verify(mockJobBuilder).submit();
  }

  @Test
  void createImputationRunStairwayError() {
    Pipeline testPipelineWithId = TestUtils.TEST_PIPELINE_1;
    testPipelineWithId.setId(1L);

    // test that when we try to create a new run but Stairway fails, we don't write the run
    // to our database (i.e. test the transaction)
    when(mockJobBuilder.flightClass(any())).thenThrow(new RuntimeException("Stairway error"));

    assertThrows(
        RuntimeException.class,
        () ->
            pipelineRunsService.createPipelineRun(
                testPipelineWithId,
                testJobId,
                testUserId,
                testDescription,
                testPipelineInputs,
                testResultUrl));

    // check that the pipeline is not written to the pipeline_runs table
    assertEquals(
        Optional.empty(), pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId));
  }
}
