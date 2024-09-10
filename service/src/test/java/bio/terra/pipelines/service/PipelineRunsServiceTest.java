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
import bio.terra.pipelines.common.utils.pagination.PageResponse;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputs;
import bio.terra.pipelines.stairway.imputation.RunImputationAzureJobFlight;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.Entity;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.net.MalformedURLException;
import java.net.URL;
import java.time.LocalDateTime;
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
  @MockBean private GcsService mockGcsService;
  @MockBean private SamService mockSamService;

  private final String testUserId = TestUtils.TEST_USER_ID_1;

  private final CommonPipelineRunStatusEnum preparingStatus = CommonPipelineRunStatusEnum.PREPARING;
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final String testDescription = TestUtils.TEST_PIPELINE_DESCRIPTION_1;
  private final String testResultUrl = TestUtils.TEST_RESULT_URL;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final String testControlWorkspaceProject = TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT;
  private final String testControlWorkspaceName = TestUtils.CONTROL_WORKSPACE_NAME;
  private final String testControlWorkspaceStorageContainerName =
      TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME;
  private final String testControlWorkspaceGoogleProject =
      TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT;

  private SimpleMeterRegistry meterRegistry;

  private PipelineRun createNewRunWithJobId(UUID jobId) {
    return createNewRunWithJobIdAndUser(jobId, testUserId);
  }

  private PipelineRun createNewRunWithJobIdAndUser(UUID jobId, String userId) {
    return new PipelineRun(
        jobId,
        userId,
        testPipelineId,
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
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
            TestUtils.TEST_PIPELINE_1.getWdlMethodVersion(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceId(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceBillingProject(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceName(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceStorageContainerName(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceGoogleProject(),
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
            testControlWorkspaceProject,
            testControlWorkspaceName,
            testControlWorkspaceStorageContainerName,
            testControlWorkspaceGoogleProject,
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
    assertEquals(testControlWorkspaceProject, savedRun.getWorkspaceBillingProject());
    assertEquals(testControlWorkspaceName, savedRun.getWorkspaceName());
    assertEquals(
        testControlWorkspaceStorageContainerName, savedRun.getWorkspaceStorageContainerName());
    assertEquals(testControlWorkspaceGoogleProject, savedRun.getWorkspaceGoogleProject());

    // verify info written to pipelineInputs table
    Optional<PipelineInput> pipelineInput = pipelineInputsRepository.findById(savedRun.getId());
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
    // missing workspace project
    Pipeline testPipelineWithIdMissingWorkspaceProject = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceProject.setWorkspaceBillingProject(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithIdMissingWorkspaceProject,
                testJobId,
                testUserId,
                testPipelineInputs));

    // missing workspace name
    Pipeline testPipelineWithIdMissingWorkspaceName = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceName.setWorkspaceName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithIdMissingWorkspaceName, testJobId, testUserId, testPipelineInputs));

    // missing workspace storage container url
    Pipeline testPipelineWithIdMissingWorkspaceStorageContainerUrl = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceStorageContainerUrl.setWorkspaceStorageContainerName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithIdMissingWorkspaceStorageContainerUrl,
                testJobId,
                testUserId,
                testPipelineInputs));
  }

  @Test
  void preparePipelineRunAlreadyExistsSameUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        testUserId,
        testPipelineId,
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
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
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
        testPipelineInputs);

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithId, testJobId, testUserId, testPipelineInputs));
  }

  @Test
  void preparePipelineRunImputation() throws MalformedURLException {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));

    Counter counter = meterRegistry.find("teaspoons.pipeline.prepare.count").counter();
    assertNull(counter);

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generatePutObjectSignedUrl(any(), any(), any())).thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineRunsService.preparePipelineRun(
            testPipelineWithId, testJobId, testUserId, userPipelineInputs);

    assertEquals(userPipelineInputs.size(), formattedPipelineFileInputs.size());
    assertEquals(
        fakeUrl.toString(), formattedPipelineFileInputs.get(fileInputKeyName).get("signedUrl"));
    assertEquals(
        "curl -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s'"
            .formatted(fileInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(fileInputKeyName).get("curlCommand"));

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
    Optional<PipelineInput> pipelineInput = pipelineInputsRepository.findById(savedRun.getId());
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
    // missing workspace project
    Pipeline testPipelineWithIdMissingWorkspaceProject = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceProject.setWorkspaceBillingProject(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithIdMissingWorkspaceProject,
                testJobId,
                testUserId,
                testDescription,
                testResultUrl));

    // missing workspace name
    Pipeline testPipelineWithIdMissingWorkspaceName = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceName.setWorkspaceName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithIdMissingWorkspaceName,
                testJobId,
                testUserId,
                testDescription,
                testResultUrl));

    // missing workspace storage container url
    Pipeline testPipelineWithIdMissingWorkspaceStorageContainerUrl = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceStorageContainerUrl.setWorkspaceStorageContainerName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithIdMissingWorkspaceStorageContainerUrl,
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
            testControlWorkspaceProject,
            testControlWorkspaceName,
            testControlWorkspaceStorageContainerName,
            testControlWorkspaceGoogleProject,
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
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
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
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
        testPipelineInputs);

    // override this mock to ensure the correct flight class is being requested
    when(mockJobBuilder.flightClass(RunImputationAzureJobFlight.class)).thenReturn(mockJobBuilder);

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
  void extractPipelineOutputsFromEntity() {
    // test that the method correctly extracts the outputs from the entity
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(
        Map.of("output_name", "gs://bucket/file1", "testNonOutputKey", "doesn't matter"));

    Map<String, String> extractedOutputs =
        pipelineRunsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);

    assertEquals(1, extractedOutputs.size());
    // the meethod should also have converted the wdlVariableName key to the camelCase outputName
    // key
    assertEquals("gs://bucket/file1", extractedOutputs.get("outputName"));
  }

  @Test
  void extractPipelineOutputsFromEntityMissingOutput() {
    // test that the method correctly throws an error if an output is missing
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(Map.of("testNonOutputKey", "doesn't matter"));

    assertThrows(
        InternalServerErrorException.class,
        () -> pipelineRunsService.extractPipelineOutputsFromEntity(outputDefinitions, entity));
  }

  @Test
  void extractPipelineOutputsFromEntityEmptyOutput() {
    // test that the method correctly throws an error if an output is empty
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(Map.of("outputName", ""));

    assertThrows(
        InternalServerErrorException.class,
        () -> pipelineRunsService.extractPipelineOutputsFromEntity(outputDefinitions, entity));
  }

  @Test
  void formatPipelineRunOutputs() throws MalformedURLException {
    PipelineRun pipelineRun = createNewRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(pipelineRun.getId());
    pipelineOutput.setOutputs(
        pipelineRunsService.pipelineRunOutputsAsString(TestUtils.TEST_PIPELINE_OUTPUTS));
    pipelineOutputsRepository.save(pipelineOutput);

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");
    // mock GCS service
    when(mockGcsService.generateGetObjectSignedUrl(any(), any(), any())).thenReturn(fakeUrl);

    ApiPipelineRunOutputs apiPipelineRunOutputs =
        pipelineRunsService.formatPipelineRunOutputs(pipelineRun);

    assertEquals(fakeUrl.toString(), apiPipelineRunOutputs.get("testFileOutputKey"));
  }

  @Test
  void markPipelineRunSuccess() {
    PipelineRun pipelineRun = createNewRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    PipelineRun updatedPipelineRun =
        pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(
            testJobId, testUserId, TestUtils.TEST_PIPELINE_OUTPUTS);
    assertTrue(updatedPipelineRun.getIsSuccess());

    Map<String, String> extractedOutputs =
        pipelineRunsService.pipelineRunOutputsAsMap(
            pipelineOutputsRepository
                .findPipelineOutputsByJobId(updatedPipelineRun.getId())
                .getOutputs());

    for (Map.Entry<String, String> entry : TestUtils.TEST_PIPELINE_OUTPUTS.entrySet()) {
      assertEquals(entry.getValue(), extractedOutputs.get(entry.getKey()));
    }
  }

  @Test
  void findPipelineRunsPaginatedNoResults() {
    PageResponse<List<PipelineRun>> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(10, null, "userIdDoesntHaveRecords");

    assertTrue(pageResults.content().isEmpty());
    assertNull(pageResults.nextPageCursor());
    assertNull(pageResults.previousPageCursor());
  }

  @Test
  void findPipelineRunsPaginatedNoOtherPages() {
    PageResponse<List<PipelineRun>> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(10, null, testUserId);

    assertEquals(1, pageResults.content().size());
    assertNull(pageResults.nextPageCursor());
    assertNull(pageResults.previousPageCursor());
  }

  @Test
  void findPageResultsResultsUseNextPage() {
    // add 3 new jobs so 4 total exist in database
    PipelineRun pipelineRun = createNewRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);
    pipelineRun = createNewRunWithJobId(UUID.randomUUID());
    pipelineRunsRepository.save(pipelineRun);
    pipelineRun = createNewRunWithJobId(UUID.randomUUID());
    pipelineRunsRepository.save(pipelineRun);

    // query for first (default) page with page size 2 so there is a next page token that exists
    PageResponse<List<PipelineRun>> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(2, null, testUserId);

    assertEquals(2, pageResults.content().size());
    assertNotNull(pageResults.nextPageCursor());
    assertNull(pageResults.previousPageCursor());
    LocalDateTime firstResultTime = pageResults.content().get(0).getCreated();

    // now query for next page
    pageResults =
        pipelineRunsService.findPipelineRunsPaginated(2, pageResults.nextPageCursor(), testUserId);

    assertEquals(2, pageResults.content().size());
    assertNull(pageResults.nextPageCursor());
    assertNotNull(pageResults.previousPageCursor());

    LocalDateTime thirdResultTime = pageResults.content().get(0).getCreated();
    // test that results are coming with most recent first
    assertTrue(firstResultTime.isAfter(thirdResultTime));
  }
}
