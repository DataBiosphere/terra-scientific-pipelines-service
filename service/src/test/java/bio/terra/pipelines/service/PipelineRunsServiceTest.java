package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

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

  private final String testStatus = "RUNNING";
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private PipelineRun createTestJobWithJobId(UUID jobId, String runStatus) {
    return createTestJobWithJobIdAndUser(jobId, testUserId, runStatus);
  }

  private PipelineRun createTestJobWithJobIdAndUser(UUID jobId, String userId, String runStatus) {
    return new PipelineRun(jobId, userId, testPipelineId, runStatus);
  }

  @Test
  void writeJobToDbOk() {
    List<PipelineRun> jobsDefault = pipelineRunsRepository.findAllByUserId(testUserId);

    // test data migration inserts one row by default
    assertEquals(1, jobsDefault.size());

    UUID savedUUID =
        pipelineRunsService.writePipelineRunToDb(
            testJobId, testUserId, testPipelineId, testPipelineInputs);

    List<PipelineRun> jobsAfterSave = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(2, jobsAfterSave.size());

    // verify info written to the jobs table
    PipelineRun savedJob =
        pipelineRunsRepository.findJobByJobIdAndUserId(savedUUID, testUserId).orElseThrow();
    assertEquals(testJobId, savedJob.getJobId());
    assertEquals(testPipelineId, savedJob.getPipelineId());
    assertEquals(testUserId, savedJob.getUserId());

    // verify info written to pipelineInputs table
    Optional<PipelineInput> pipelineInput = pipelineInputsRepository.findById(savedJob.getId());
    assertTrue(pipelineInput.isPresent());
    assertEquals("{\"first_key\":\"first_value\"}", pipelineInput.get().getInputs());
  }

  @Test
  void writeJobToDbDuplicateJob() {
    // try to save a job with the same job id two times, the second time it should throw duplicate
    // exception error
    PipelineRun newPipelineRun = createTestJobWithJobId(testJobId, testStatus);

    PipelineRun savedJobFirst =
        pipelineRunsService.writePipelineRunToDbThrowsDuplicateException(newPipelineRun);
    assertNotNull(savedJobFirst);

    PipelineRun newPipelineRunSameId = createTestJobWithJobId(testJobId, testStatus);
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
    PipelineRun newPipelineRun = createTestJobWithJobId(testJobId, testStatus);

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
        createTestJobWithJobIdAndUser(UUID.randomUUID(), testUserId2, testStatus);
    pipelineRunsRepository.save(newPipelineRun);

    // Verify that the old userid still show only 1 record
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId);
    assertEquals(1, pipelineRuns.size());

    // Verify the new user's id shows a single job as well
    pipelineRuns = pipelineRunsRepository.findAllByUserId(testUserId2);
    assertEquals(1, pipelineRuns.size());
  }
}
