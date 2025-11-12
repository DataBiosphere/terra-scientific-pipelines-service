package bio.terra.pipelines.service;

import static bio.terra.pipelines.service.PipelineRunsService.ALLOWED_SORT_PROPERTIES;
import static bio.terra.pipelines.testutils.TestUtils.createNewPipelineRunWithJobId;
import static bio.terra.pipelines.testutils.TestUtils.createNewPipelineRunWithJobIdAndUser;
import static bio.terra.pipelines.testutils.TestUtils.createTestPipelineWithId;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.pagination.PageResponse;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.flights.imputation.v20251002.RunImputationGcpJobFlight;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.Flight;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.net.MalformedURLException;
import java.net.URL;
import java.time.Instant;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.domain.Page;
import org.springframework.data.domain.Pageable;
import org.springframework.data.domain.Sort;
import org.springframework.data.jpa.domain.Specification;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class PipelineRunsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelineRunsService pipelineRunsService;
  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @Autowired PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired PipelineOutputsRepository pipelineOutputsRepository;

  // mock Stairway and other services
  @MockitoBean private JobService mockJobService;
  @MockitoBean private JobBuilder mockJobBuilder;
  @MockitoBean private GcsService mockGcsService;
  @MockitoBean private SamService mockSamService;

  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final String testToolVersion = TestUtils.TEST_TOOL_VERSION_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final String testControlWorkspaceProject = TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT;
  private final String testControlWorkspaceName = TestUtils.CONTROL_WORKSPACE_NAME;
  private final String testControlWorkspaceStorageContainerName =
      TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME;
  private final String testControlWorkspaceGoogleProject =
      TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT;
  private final String testUserDescription = TestUtils.TEST_USER_PROVIDED_DESCRIPTION;
  private final Integer testQuotaConsumed = 10;
  private final Integer testRawQuotaConsumed = 5;

  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any(UUID.class))).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(Flight.class)).thenReturn(mockJobBuilder);
    when(mockJobBuilder.addParameter(anyString(), any())).thenReturn(mockJobBuilder);
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
            testToolVersion,
            testControlWorkspaceProject,
            testControlWorkspaceName,
            testControlWorkspaceStorageContainerName,
            testControlWorkspaceGoogleProject,
            testPipelineInputs,
            testUserDescription);

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
    assertEquals(testUserDescription, savedRun.getDescription());
    assertEquals(CommonPipelineRunStatusEnum.PREPARING, savedRun.getStatus());
    assertNotNull(savedRun.getCreated());
    assertNotNull(savedRun.getUpdated());
    assertEquals(testControlWorkspaceProject, savedRun.getWorkspaceBillingProject());
    assertEquals(testControlWorkspaceName, savedRun.getWorkspaceName());
    assertEquals(
        testControlWorkspaceStorageContainerName, savedRun.getWorkspaceStorageContainerName());
    assertEquals(testControlWorkspaceGoogleProject, savedRun.getWorkspaceGoogleProject());

    // verify info written to pipeline_inputs table
    Map<String, Object> pipelineInputs =
        pipelineInputsOutputsService.retrievePipelineInputs(savedRun);
    assertNotNull(pipelineInputs);
    assertEquals(
        "{\"first_key\":\"first_value\"}",
        pipelineInputsOutputsService.mapToString(pipelineInputs));
  }

  @Test
  void writeRunToDbDuplicateRun() {
    // try to save a run with the same job id two times, the second time it should throw duplicate
    // exception error
    PipelineRun newPipelineRun = createNewPipelineRunWithJobId(testJobId);

    PipelineRun savedJobFirst =
        pipelineRunsService.writePipelineRunToDbThrowsDuplicateException(newPipelineRun);
    assertNotNull(savedJobFirst);

    PipelineRun newPipelineRunSameId = createNewPipelineRunWithJobId(testJobId);
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
    PipelineRun newPipelineRun = createNewPipelineRunWithJobId(testJobId);

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
    PipelineRun newPipelineRun =
        createNewPipelineRunWithJobIdAndUser(UUID.randomUUID(), testUserId2);
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
                testPipelineInputs,
                testUserDescription));

    // missing workspace name
    Pipeline testPipelineWithIdMissingWorkspaceName = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceName.setWorkspaceName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithIdMissingWorkspaceName,
                testJobId,
                testUserId,
                testPipelineInputs,
                testUserDescription));

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
                testPipelineInputs,
                testUserDescription));
  }

  @Test
  void preparePipelineRunAlreadyExistsSameUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        testUserId,
        testPipelineId,
        testToolVersion,
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
        testPipelineInputs,
        testUserDescription);

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithId,
                testJobId,
                testUserId,
                testPipelineInputs,
                testUserDescription));
  }

  @Test
  void preparePipelineRunAlreadyExistsOtherUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        TestUtils.TEST_USER_ID_2, // different user than the caller
        testPipelineId,
        testToolVersion,
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
        testPipelineInputs,
        testUserDescription);

    assertThrows(
        BadRequestException.class,
        () ->
            pipelineRunsService.preparePipelineRun(
                testPipelineWithId,
                testJobId,
                testUserId,
                testPipelineInputs,
                testUserDescription));
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

    when(mockGcsService.generatePutObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceGoogleProject()),
            eq(testPipelineWithId.getWorkspaceStorageContainerName()),
            anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineRunsService.preparePipelineRun(
            testPipelineWithId, testJobId, testUserId, userPipelineInputs, testUserDescription);

    assertEquals(userPipelineInputs.size(), formattedPipelineFileInputs.size());
    assertEquals(
        fakeUrl.toString(), formattedPipelineFileInputs.get(fileInputKeyName).get("signedUrl"));
    assertEquals(
        "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s' | cat"
            .formatted(fileInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(fileInputKeyName).get("curlCommand"));

    // check db for the pipeline run
    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId).orElseThrow();

    assertEquals(testJobId, writtenPipelineRun.getJobId());
    assertEquals(testUserId, writtenPipelineRun.getUserId());
    assertEquals(testUserDescription, writtenPipelineRun.getDescription());
    assertEquals(testPipelineWithId.getId(), writtenPipelineRun.getPipelineId());
    assertEquals(testPipelineWithId.getToolVersion(), writtenPipelineRun.getToolVersion());
    assertEquals(testUserDescription, writtenPipelineRun.getDescription());
    assertNotNull(writtenPipelineRun.getCreated());
    assertNotNull(writtenPipelineRun.getUpdated());

    // verify info written to pipeline_inputs table
    Map<String, Object> pipelineInputs =
        pipelineInputsOutputsService.retrievePipelineInputs(writtenPipelineRun);
    assertNotNull(pipelineInputs);
    assertEquals(
        "{\"%s\":\"%s\"}".formatted(fileInputKeyName, fileInputValue),
        pipelineInputsOutputsService.mapToString(pipelineInputs));

    // verify the pipeline prepareRun counter was incremented
    counter = meterRegistry.find("teaspoons.pipeline.prepareRun.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }

  @Test
  void preparePipelineRunImputationAllowNullDescription() throws MalformedURLException {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));

    Counter counter = meterRegistry.find("teaspoons.pipeline.prepare.count").counter();
    assertNull(counter);

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generatePutObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceGoogleProject()),
            eq(testPipelineWithId.getWorkspaceStorageContainerName()),
            anyString()))
        .thenReturn(fakeUrl);

    pipelineRunsService.preparePipelineRun(
        testPipelineWithId, testJobId, testUserId, userPipelineInputs, null);

    // check db for the pipeline run
    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId).orElseThrow();

    assertEquals(testJobId, writtenPipelineRun.getJobId());
    assertEquals(testUserId, writtenPipelineRun.getUserId());
    assertNull(writtenPipelineRun.getDescription());
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
                testPipelineWithIdMissingWorkspaceProject, testJobId, testUserId));

    // missing workspace name
    Pipeline testPipelineWithIdMissingWorkspaceName = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceName.setWorkspaceName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithIdMissingWorkspaceName, testJobId, testUserId));

    // missing workspace storage container url
    Pipeline testPipelineWithIdMissingWorkspaceStorageContainerUrl = createTestPipelineWithId();
    testPipelineWithIdMissingWorkspaceStorageContainerUrl.setWorkspaceStorageContainerName(null);

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineRunsService.startPipelineRun(
                testPipelineWithIdMissingWorkspaceStorageContainerUrl, testJobId, testUserId));
  }

  @Test
  void startPipelineRunNoPreparedPipelineRun() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    assertThrows(
        BadRequestException.class,
        () -> pipelineRunsService.startPipelineRun(testPipelineWithId, testJobId, testUserId));
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
            testToolVersion,
            testControlWorkspaceProject,
            testControlWorkspaceName,
            testControlWorkspaceStorageContainerName,
            testControlWorkspaceGoogleProject,
            CommonPipelineRunStatusEnum.RUNNING,
            testUserDescription));

    assertThrows(
        BadRequestException.class,
        () -> pipelineRunsService.startPipelineRun(testPipelineWithId, testJobId, testUserId));
  }

  @Test
  void startPipelineRunWrongUser() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        TestUtils.TEST_USER_ID_2, // different user than the caller
        testPipelineId,
        testToolVersion,
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
        testPipelineInputs,
        testUserDescription);

    assertThrows(
        BadRequestException.class,
        () -> pipelineRunsService.startPipelineRun(testPipelineWithId, testJobId, testUserId));
  }

  @Test
  void startPipelineRunImputation() {
    Pipeline testPipelineWithId = createTestPipelineWithId();

    // write a prepared pipeline run to the db
    pipelineRunsService.writeNewPipelineRunToDb(
        testJobId,
        testUserId,
        testPipelineId,
        testToolVersion,
        testControlWorkspaceProject,
        testControlWorkspaceName,
        testControlWorkspaceStorageContainerName,
        testControlWorkspaceGoogleProject,
        testPipelineInputs,
        testUserDescription);

    // override this mock to ensure the correct flight class is being requested
    when(mockJobBuilder.flightClass(RunImputationGcpJobFlight.class)).thenReturn(mockJobBuilder);

    PipelineRun returnedPipelineRun =
        pipelineRunsService.startPipelineRun(testPipelineWithId, testJobId, testUserId);

    assertEquals(testJobId, returnedPipelineRun.getJobId());
    assertEquals(testUserId, returnedPipelineRun.getUserId());
    assertEquals(testUserDescription, returnedPipelineRun.getDescription());
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
    when(mockJobBuilder.flightClass(RunImputationGcpJobFlight.class))
        .thenThrow(new RuntimeException("Stairway error"));

    assertThrows(
        RuntimeException.class,
        () -> pipelineRunsService.startPipelineRun(testPipelineWithId, testJobId, testUserId));

    // check that the pipeline is not persisted to the pipeline_runs table
    assertEquals(
        Optional.empty(), pipelineRunsRepository.findByJobIdAndUserId(testJobId, testUserId));
  }

  @Test
  void markPipelineRunSuccess() {
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    PipelineRun updatedPipelineRun =
        pipelineRunsService.markPipelineRunSuccessAndWriteOutputs(
            testJobId,
            testUserId,
            TestUtils.TEST_PIPELINE_OUTPUTS,
            testQuotaConsumed,
            testRawQuotaConsumed);
    assertTrue(updatedPipelineRun.getStatus().isSuccess());
    assertEquals(CommonPipelineRunStatusEnum.SUCCEEDED, updatedPipelineRun.getStatus());
    assertEquals(testQuotaConsumed, updatedPipelineRun.getQuotaConsumed());
    assertEquals(testRawQuotaConsumed, updatedPipelineRun.getRawQuotaConsumed());

    PipelineOutput pipelineOutput =
        pipelineOutputsRepository.findPipelineOutputsByJobId(pipelineRun.getId());
    Map<String, Object> extractedOutput =
        pipelineInputsOutputsService.stringToMap(pipelineOutput.getOutputs());

    for (Map.Entry<String, String> entry : TestUtils.TEST_PIPELINE_OUTPUTS.entrySet()) {
      assertEquals(entry.getValue(), extractedOutput.get(entry.getKey()));
    }
  }

  @Nested
  @DisplayName("findPipelineRunsPaginated v2 tests")
  class FindPipelineRunsPaginatedV2Tests {
    @Test
    void findPipelineRunsPaginatedNoResultsV2() {
      Page<PipelineRun> pageResults =
          pipelineRunsService.findPipelineRunsPaginated(
              1, 10, "created", "DESC", "userIdDoesntHaveRecords");

      assertTrue(pageResults.stream().toList().isEmpty());
    }

    @Test
    void findPageResultsResultsUseNextPageV2() {
      // add 3 new jobs so 4 total exist in database
      PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
      pipelineRunsRepository.save(pipelineRun);
      pipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
      pipelineRunsRepository.save(pipelineRun);
      pipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
      pipelineRunsRepository.save(pipelineRun);

      // query for first (default) page with page size 2 so there is a next page token that exists
      Page<PipelineRun> pageResults =
          pipelineRunsService.findPipelineRunsPaginated(0, 2, null, null, testUserId);

      assertEquals(2, pageResults.stream().toList().size());
      Instant firstResultTime = pageResults.stream().findFirst().get().getCreated();

      // now query for next page
      pageResults = pipelineRunsService.findPipelineRunsPaginated(1, 2, null, null, testUserId);

      assertEquals(2, pageResults.stream().toList().size());

      Instant thirdResultTime = pageResults.stream().findFirst().get().getCreated();
      // test that results are coming with most recent first
      assertTrue(firstResultTime.isAfter(thirdResultTime));
    }

    @Test
    void findPipelineRunsPaginatedAndSortedV2() {
      PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
      pipelineRunsRepository.save(pipelineRun);
      pipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
      pipelineRunsRepository.save(pipelineRun);
      pipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
      pipelineRunsRepository.save(pipelineRun);

      Page<PipelineRun> pageResultsAsc =
          pipelineRunsService.findPipelineRunsPaginated(0, 10, "created", "ASC", testUserId);

      Page<PipelineRun> pageResultsDesc =
          pipelineRunsService.findPipelineRunsPaginated(0, 10, "created", "DESC", testUserId);

      assertEquals(4, pageResultsAsc.stream().toList().size());
      assertEquals(4, pageResultsDesc.stream().toList().size());

      // check that the ascending and descending results are the reverse of each other
      List<PipelineRun> ascList = pageResultsAsc.stream().toList();
      List<PipelineRun> descList = pageResultsDesc.stream().toList();
      for (int i = 0; i < ascList.size(); i++) {
        assertEquals(ascList.get(i).getId(), descList.get(descList.size() - 1 - i).getId());
      }
    }

    @Test
    void findPipelineRunsPaginatedDefaultSortDirectionV2() {
      PipelineRunsRepository mockPipelineRunsRepository =
          org.mockito.Mockito.mock(PipelineRunsRepository.class);
      ArgumentCaptor<Pageable> pageableCaptor = ArgumentCaptor.forClass(Pageable.class);
      ArgumentCaptor<Specification<PipelineRun>> specCaptor =
          ArgumentCaptor.forClass(Specification.class);

      PipelineRunsService mockPipelineRunsService =
          new PipelineRunsService(
              mockJobService, pipelineInputsOutputsService, mockPipelineRunsRepository, null, null);

      // query with null sort params, should default to created DESC
      mockPipelineRunsService.findPipelineRunsPaginated(0, 10, "created", null, testUserId);

      verify(mockPipelineRunsRepository).findAll(specCaptor.capture(), pageableCaptor.capture());

      // Assert the sort direction is 'DESC'
      Pageable capturedPageable = pageableCaptor.getValue();
      assertTrue(capturedPageable.getSort().isSorted());
      assertEquals(
          Sort.Direction.DESC, capturedPageable.getSort().getOrderFor("created").getDirection());
    }

    @Test
    void findPipelineRunsPaginatedDefaultSortPropertyV2() {
      PipelineRunsRepository mockPipelineRunsRepository =
          org.mockito.Mockito.mock(PipelineRunsRepository.class);
      ArgumentCaptor<Pageable> pageableCaptor = ArgumentCaptor.forClass(Pageable.class);
      ArgumentCaptor<Specification<PipelineRun>> specCaptor =
          ArgumentCaptor.forClass(Specification.class);

      PipelineRunsService mockPipelineRunsService =
          new PipelineRunsService(
              mockJobService, pipelineInputsOutputsService, mockPipelineRunsRepository, null, null);

      // query with null sort property, should default to created
      mockPipelineRunsService.findPipelineRunsPaginated(0, 10, null, "DESC", testUserId);

      verify(mockPipelineRunsRepository).findAll(specCaptor.capture(), pageableCaptor.capture());

      // Assert the sort property is 'created'
      Pageable capturedPageable = pageableCaptor.getValue();
      assertTrue(capturedPageable.getSort().isSorted());
      assertEquals("created", capturedPageable.getSort().getOrderFor("created").getProperty());
    }

    @Test
    void findPipelineRunsPaginatedInvalidSortPropertyV2() {
      assertThrows(
          BadRequestException.class,
          () -> pipelineRunsService.findPipelineRunsPaginated(0, 10, "INVALID", "ASC", testUserId));
    }

    @Test
    void findPipelineRunsPaginatedInvalidSortDirectionV2() {
      assertThrows(
          BadRequestException.class,
          () ->
              pipelineRunsService.findPipelineRunsPaginated(
                  0, 10, "created", "INVALID", testUserId));
    }
  }

  @Nested
  @Deprecated
  @DisplayName("findPipelineRunsPaginated v1 tests")
  class FindPipelineRunsPaginatedV1Tests {

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
      PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
      pipelineRunsRepository.save(pipelineRun);
      pipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
      pipelineRunsRepository.save(pipelineRun);
      pipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
      pipelineRunsRepository.save(pipelineRun);

      // query for first (default) page with page size 2 so there is a next page token that exists
      PageResponse<List<PipelineRun>> pageResults =
          pipelineRunsService.findPipelineRunsPaginated(2, null, testUserId);

      assertEquals(2, pageResults.content().size());
      assertNotNull(pageResults.nextPageCursor());
      assertNull(pageResults.previousPageCursor());
      Instant firstResultTime = pageResults.content().get(0).getCreated();

      // now query for next page
      pageResults =
          pipelineRunsService.findPipelineRunsPaginated(
              2, pageResults.nextPageCursor(), testUserId);

      assertEquals(2, pageResults.content().size());
      assertNull(pageResults.nextPageCursor());
      assertNotNull(pageResults.previousPageCursor());

      Instant thirdResultTime = pageResults.content().get(0).getCreated();
      // test that results are coming with most recent first
      assertTrue(firstResultTime.isAfter(thirdResultTime));
    }
  }

  @Test
  void getPipelineRunCountForUser() {
    // A test row already exists for this user. Confirm the initial count is 1.
    long initialResultsCount = pipelineRunsService.getPipelineRunCount(testUserId);
    assertEquals(1, initialResultsCount);

    // Add a new pipeline run for the same user and confirm the count is now 2.
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    long totalResultsCount = pipelineRunsService.getPipelineRunCount(testUserId);
    assertEquals(2, totalResultsCount);
  }

  @Test
  void getFilteredPipelineRunCountForUser() {
    // A test row already exists for this user with status RUNNING. Confirm the initial count.
    Map<String, String> filtersForStatus = new HashMap<>();
    filtersForStatus.put("status", "RUNNING");
    long initialResultsCount =
        pipelineRunsService.getFilteredPipelineRunCount(testUserId, filtersForStatus);
    assertEquals(1, initialResultsCount);

    // Add a new pipeline run for the same user with status SUCCEEDED
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRunsRepository.save(pipelineRun);

    long preparingResultsCount =
        pipelineRunsService.getFilteredPipelineRunCount(testUserId, filtersForStatus);
    assertEquals(1, preparingResultsCount);

    filtersForStatus.put("status", "SUCCEEDED");
    long succeededResultsCount =
        pipelineRunsService.getFilteredPipelineRunCount(testUserId, filtersForStatus);
    assertEquals(1, succeededResultsCount);

    // Confirm total count without filters
    long totalResultsCount =
        pipelineRunsService.getFilteredPipelineRunCount(testUserId, new HashMap<>());
    assertEquals(2, totalResultsCount);
  }

  @Test
  void validateAllowedSortPropertiesExist() {
    // This test enforces that the allowed sort properties are kept
    // up to date if the PipelineRun class changes.
    for (String property : ALLOWED_SORT_PROPERTIES) {
      assertDoesNotThrow(
          () -> PipelineRun.class.getDeclaredField(property),
          "Sort property '%s' was expected to exist on PipelineRun class".formatted(property));
    }
  }
}
