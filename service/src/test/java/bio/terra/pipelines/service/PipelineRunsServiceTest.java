package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputs;
import bio.terra.pipelines.db.entities.PipelineOutputs;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutput;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlight;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.HashMap;
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
  @Autowired PipelineOutputsRepository pipelineOutputsRepository;

  // mock Stairway and other services
  @MockBean private JobService mockJobService;
  @MockBean private JobBuilder mockJobBuilder;
  @MockBean private WorkspaceManagerService mockWorkspaceManagerService;
  @MockBean private SamService mockSamService;

  private final String testUserId = TestUtils.TEST_USER_ID_1;

  private final CommonPipelineRunStatusEnum preparingStatus = CommonPipelineRunStatusEnum.PREPARING;
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final String testDescription = TestUtils.TEST_PIPELINE_DESCRIPTION_1;
  private final String testResultUrl = TestUtils.TEST_RESULT_URL;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final UUID testControlWorkspaceId = TestUtils.CONTROL_WORKSPACE_ID;
  private final String testControlWorkspaceStorageContainerUrl =
      TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_URL;

  private SimpleMeterRegistry meterRegistry;

  private PipelineRun createNewRunWithJobId(UUID jobId) {
    return createNewRunWithJobIdAndUser(jobId, testUserId);
  }

  private PipelineRun createNewRunWithJobIdAndUser(UUID jobId, String userId) {
    return new PipelineRun(
        jobId,
        userId,
        testPipelineId,
        testControlWorkspaceId,
        testControlWorkspaceStorageContainerUrl,
        preparingStatus);
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
            TestUtils.TEST_PIPELINE_1.getPipelineInputDefinitions(),
            TestUtils.TEST_PIPELINE_1.getPipelineOutputDefinitions());
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

    when(mockSamService.getTeaspoonsServiceAccountToken()).thenReturn("teaspoonsSaToken");

    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.add(meterRegistry);
  }

  @Test
  void writeNewRunToDbOk() {
    List<PipelineRun> runsDefault = pipelineRunsRepository.findAllByUserId(testUserId);

    // test data migration inserts one row by default
    assertEquals(1, runsDefault.size());

    PipelineRun pipelineRun =
        pipelineRunsService.writeNewPipelineRunToDb(
            testJobId,
            testUserId,
            testPipelineId,
            testControlWorkspaceId,
            testControlWorkspaceStorageContainerUrl,
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
    assertNull(savedRun.getDescription());
    assertNull(savedRun.getResultUrl());
    assertEquals(CommonPipelineRunStatusEnum.PREPARING, savedRun.getStatus());
    assertNotNull(savedRun.getCreated());
    assertNotNull(savedRun.getUpdated());

    // verify info written to pipelineInputs table
    Optional<PipelineInputs> pipelineInput = pipelineInputsRepository.findById(savedRun.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals("{\"first_key\":\"first_value\"}", pipelineInput.get().getInputs());
  }

  @Test
  void writeRunToDbDuplicateRun() {
    // try to save a run with the same job id two times, the second time it should throw duplicate
    // exception error
    PipelineRun newPipelineRun = createNewRunWithJobId(testJobId);

    PipelineRun savedJobFirst =
        pipelineRunsService.writePipelineRunToDbThrowsDuplicateException(newPipelineRun);
    assertNotNull(savedJobFirst);

    PipelineRun newPipelineRunSameId = createNewRunWithJobId(testJobId);
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
    PipelineRun newPipelineRun = createNewRunWithJobId(testJobId);

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
    PipelineRun newPipelineRun = createNewRunWithJobIdAndUser(UUID.randomUUID(), testUserId2);
    pipelineRunsRepository.save(newPipelineRun);

    // Verify that the old userid still show only 1 record
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(1, pipelineRuns.size());

    // Verify the new user's id shows a single job as well
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId2);
    assertEquals(1, pipelineRuns.size());
  }

  @Test
  void preparePipelineRunNoWorkspaceSetUp() {
    Pipeline testPipelineWithIdMissingWorkspaceId = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceId.setWorkspaceId(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithIdMissingWorkspaceId, testJobId, testUserId, testPipelineInputs));
  }

  @Test
  void preparePipelineRunAlreadyExistsSameUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        testUserId,
        testPipelineId,
        testControlWorkspaceId,
        testControlWorkspaceStorageContainerUrl,
        testPipelineInputs);

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithId, testJobId, testUserId, testPipelineInputs));
  }

  @Test
  void preparePipelineRunAlreadyExistsOtherUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        TestUtils.TEST_USER_ID_2, // different user than the caller
        testPipelineId,
        testControlWorkspaceId,
        testControlWorkspaceStorageContainerUrl,
        testPipelineInputs);

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithId, testJobId, testUserId, testPipelineInputs));
  }

  @Test
  void preparePipelineRunImputation() {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));
    String sasUrlValue =
        "https://lz123.blob.core.windows.net/sc-%s/fake/file.vcf.gz?somestuff"
            .formatted(testPipelineWithId.getWorkspaceId());

    Counter counter = meterRegistry.find("teaspoons.pipeline.prepare.count").counter();
    assertNull(counter);

    // mocks
    when(mockWorkspaceManagerService.getWorkspaceStorageResourceId(any(), any()))
        .thenReturn(UUID.randomUUID());
    when(mockWorkspaceManagerService.getWriteSasUrlForBlob(any(), any(), any(), any()))
        .thenReturn(sasUrlValue);

    Map<String, Map<String, String>> fileUploadsMap =
        pipelineRunsService.preparePipelineRun(
            testPipelineWithId, testJobId, testUserId, userPipelineInputs);

    // test input definitions have one user-provided file input
    assertEquals(1, fileUploadsMap.size());
    assertEquals(sasUrlValue, fileUploadsMap.get(fileInputKeyName).get("sasUrl"));
    assertEquals(
        "azcopy copy %s %s".formatted(fileInputValue, sasUrlValue),
        fileUploadsMap.get(fileInputKeyName).get("azcopyCommand"));

    // check db for the pipeline run
    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId).orElseThrow();

    assertEquals(testJobId, writtenPipelineRun.getJobId());
    assertEquals(testUserId, writtenPipelineRun.getUserId());
    assertNull(writtenPipelineRun.getDescription());
    assertNull(writtenPipelineRun.getResultUrl());
    assertEquals(testPipelineWithId.getId(), writtenPipelineRun.getPipelineId());
    assertNotNull(writtenPipelineRun.getCreated());
    assertNotNull(writtenPipelineRun.getUpdated());

    // verify info written to pipeline_runs table
    PipelineRun savedRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(writtenPipelineRun.getJobId(), testUserId)
            .orElseThrow();
    assertEquals(testJobId, savedRun.getJobId());
    assertNotNull(savedRun.getCreated());
    assertNotNull(savedRun.getUpdated());

    // verify info written to pipeline_inputs table
    Optional<PipelineInputs> pipelineInput = pipelineInputsRepository.findById(savedRun.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals(
        "{\"%s\":\"%s\"}".formatted(fileInputKeyName, fileInputValue),
        pipelineInput.get().getInputs());

    // verify the pipeline prepareRun counter was incremented
    counter = meterRegistry.find("teaspoons.pipeline.prepareRun.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }

  @Test
  void startPipelineRunNoWorkspaceSetUp() {
    Pipeline testPipelineWithIdMissingWorkspaceId = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceId.setWorkspaceId(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithIdMissingWorkspaceId,
                testJobId,
                testUserId,
                testDescription,
                testResultUrl));
  }

  @Test
  void startPipelineRunNoPreparedPipelineRun() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithId, testJobId, testUserId, testDescription, testResultUrl));
  }

  @Test
  void startPipelineRunAlreadyRunning() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a RUNNING pipelineRun to the db
    pipelineRunsRepository.save(
        new PipelineRun(
            testJobId,
            testUserId,
            testPipelineId,
            testControlWorkspaceId,
            testControlWorkspaceStorageContainerUrl,
            CommonPipelineRunStatusEnum.RUNNING));

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithId, testJobId, testUserId, testDescription, testResultUrl));
  }

  @Test
  void startPipelineRunWrongUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        TestUtils.TEST_USER_ID_2, // different user than the caller
        testPipelineId,
        testControlWorkspaceId,
        testControlWorkspaceStorageContainerUrl,
        testPipelineInputs);

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithId, testJobId, testUserId, testDescription, testResultUrl));
  }

  @Test
  void startPipelineRunImputation() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        testUserId,
        testPipelineId,
        testControlWorkspaceId,
        testControlWorkspaceStorageContainerUrl,
        testPipelineInputs);

    // override this mock to ensure the correct flight class is being requested
    when(mockJobBuilder.flightClass(RunImputationJobFlight.class)).thenReturn(mockJobBuilder);

    PipelineRun returnedPipelineRun =
        pipelineRunsService.startPipelineRun(
            testPipelineWithId, testJobId, testUserId, testDescription, testResultUrl);

    assertEquals(testJobId, returnedPipelineRun.getJobId());
    assertEquals(testUserId, returnedPipelineRun.getUserId());
    assertEquals(testDescription, returnedPipelineRun.getDescription());
    assertEquals(testResultUrl, returnedPipelineRun.getResultUrl());
    assertEquals(testPipelineWithId.getId(), returnedPipelineRun.getPipelineId());
    assertNotNull(returnedPipelineRun.getCreated());
    assertNotNull(returnedPipelineRun.getUpdated());

    // verify info written to pipeline_runs table
    PipelineRun savedRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(returnedPipelineRun.getJobId(), testUserId)
            .orElseThrow();
    assertEquals(testJobId, savedRun.getJobId());
    assertEquals(CommonPipelineRunStatusEnum.RUNNING, savedRun.getStatus());
    assertNotNull(savedRun.getCreated());
    assertNotNull(savedRun.getUpdated());

    // verify submit was called
    verify(mockJobBuilder).submit();
  }

  @Test
  void createImputationRunStairwayError() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // test that when we try to create a new run but Stairway fails, we don't write the run
    // to our database (i.e. test the transaction)
    when(mockJobBuilder.flightClass(any())).thenThrow(new RuntimeException("Stairway error"));

    assertThrows(
        RuntimeException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithId, testJobId, testUserId, testDescription, testResultUrl));

    // check that the pipeline is not persisted to the pipeline_runs table
    assertEquals(
        Optional.empty(), pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId));
  }

  @Test
  void formatPipelineRunOutputs() {
    PipelineRun pipelineRun = createNewRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    PipelineOutputs pipelineOutputs = new PipelineOutputs();
    pipelineOutputs.setJobId(pipelineRun.getId());
    pipelineOutputs.setOutputs(
        pipelineRunsService.pipelineRunOutputAsString(TestUtils.TEST_PIPELINE_OUTPUTS));
    pipelineOutputsRepository.save(pipelineOutputs);

    String sasUrl = "sasUrlValue";
    // mock WorkspaceManagerService
    when(mockWorkspaceManagerService.getReadSasUrlForBlob(any(), any(), any(), any()))
        .thenReturn(sasUrl);

    ApiPipelineRunOutput apiPipelineRunOutput =
        pipelineRunsService.formatPipelineRunOutputs(pipelineRun);

    assertEquals(sasUrl, apiPipelineRunOutput.get("testFileOutputKey"));
  }

  @Test
  void markPipelineRunSuccess() {
    PipelineRun pipelineRun = createNewRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    PipelineRun updatedPipelineRun =
        pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(
            testJobId, testUserId, TestUtils.TEST_PIPELINE_OUTPUTS);
    assertTrue(updatedPipelineRun.getIsSuccess());

    Map<String, String> extractedOutput =
        pipelineRunsService.pipelineRunOutputAsMap(
            pipelineOutputsRepository
                .findPipelineOutputsByJobId(updatedPipelineRun.getId())
                .getOutputs());

    for (Map.Entry<String, String> entry : TestUtils.TEST_PIPELINE_OUTPUTS.entrySet()) {
      assertEquals(entry.getValue(), extractedOutput.get(entry.getKey()));
    }
  }
}
