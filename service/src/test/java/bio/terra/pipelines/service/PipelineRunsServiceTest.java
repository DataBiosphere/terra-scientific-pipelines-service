package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class PipelineRunsServiceTest extends BaseEmbeddedDbTest {
  @Autowired PipelineRunsService pipelineRunsService;
  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @Autowired PipelineInputsRepository pipelineInputsRepository;

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
}
