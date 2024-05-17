package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class ImputationServiceTest extends BaseEmbeddedDbTest {

  @Autowired ImputationService imputationService;
  @Autowired JobsRepository jobsRepository;
  @Autowired PipelineInputsRepository pipelineInputsRepository;

  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private Job createTestJobWithJobId(UUID jobId) {
    return createTestJobWithJobIdAndUser(jobId, testUserId);
  }

  private Job createTestJobWithJobIdAndUser(UUID jobId, String userId) {
    return new Job(jobId, userId, testPipelineId);
  }

  @Test
  void writeJobToDbOk() {
    List<Job> jobsDefault = jobsRepository.findAllByUserId(testUserId);

    // test data migration inserts one row by default
    assertEquals(1, jobsDefault.size());

    UUID savedUUID =
        imputationService.writeJobToDb(testJobId, testUserId, testPipelineId, testPipelineInputs);

    List<Job> jobsAfterSave = jobsRepository.findAllByUserId(testUserId);
    assertEquals(2, jobsAfterSave.size());

    // verify info written to the jobs table
    Job savedJob = jobsRepository.findJobByJobIdAndUserId(savedUUID, testUserId).orElseThrow();
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
    Job newJob = createTestJobWithJobId(testJobId);

    Job savedJobFirst = imputationService.writeJobToDbThrowsDuplicateException(newJob);
    assertNotNull(savedJobFirst);

    Job newJobSameId = createTestJobWithJobId(testJobId);
    assertThrows(
        DuplicateObjectException.class,
        () -> imputationService.writeJobToDbThrowsDuplicateException(newJobSameId));
  }

  @Test
  void testGetCorrectNumberOfRows() {
    // A test row should exist for this user.
    List<Job> jobs = jobsRepository.findAllByUserId(testUserId);
    assertEquals(1, jobs.size());

    // insert another row and verify that it shows up
    Job newJob = createTestJobWithJobId(testJobId);

    jobsRepository.save(newJob);
    jobs = jobsRepository.findAllByUserId(testUserId);
    assertEquals(2, jobs.size());
  }

  @Test
  void testCorrectUserIsolation() {
    // A test row should exist for this user.
    List<Job> jobs = jobsRepository.findAllByUserId(testUserId);
    assertEquals(1, jobs.size());

    // insert row for second user and verify that it shows up
    String testUserId2 = TestUtils.TEST_USER_ID_2;
    Job newJob = createTestJobWithJobIdAndUser(UUID.randomUUID(), testUserId2);
    jobsRepository.save(newJob);

    // Verify that the old userid still show only 1 record
    jobs = jobsRepository.findAllByUserId(testUserId);
    assertEquals(1, jobs.size());

    // Verify the new user's id shows a single job as well
    jobs = jobsRepository.findAllByUserId(testUserId2);
    assertEquals(1, jobs.size());
  }
}
