package bio.terra.pipelines.db.repositories;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import jakarta.persistence.EntityManager;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.transaction.annotation.Transactional;

@Transactional
class PipelineInputsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Autowired private PipelineInputsRepository pipelineInputsRepository;
  @Autowired private EntityManager entityManager;

  @Test
  // This test verifies that cascade delete is configured correctly and that
  // deleting a pipeline run also deletes all associated inputs
  void deletingPipelineRunDeletesAssociatedInputs() {
    // create and save a pipeline run
    UUID jobId = UUID.randomUUID();
    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(jobId);
    pipelineRun = pipelineRunsRepository.save(pipelineRun);

    // create and save inputs
    String inputsJson = "{\"input1\":\"gs://my-input-bucket/input1\",\"input2\":\"helloWorld\"}";
    PipelineInput input = new PipelineInput();
    input.setPipelineRunId(pipelineRun.getId());
    input.setInputs(inputsJson);
    pipelineInputsRepository.save(input);

    // verify pipeline run exists
    PipelineRun retrievedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNotNull(retrievedRun);

    // verify input exists
    PipelineInput retrievedInput =
        pipelineInputsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNotNull(retrievedInput);
    assertEquals(jobId, retrievedRun.getJobId());
    assertEquals(inputsJson, retrievedInput.getInputs());

    // delete pipeline run
    pipelineRunsRepository.delete(pipelineRun);

    // flush and clear the persistence context to ensure we are fetching fresh data from the
    // database
    entityManager.flush();
    entityManager.clear();

    // verify pipeline run is deleted
    PipelineRun deletedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNull(deletedRun);

    // verify cascade delete - input should also be deleted
    PipelineInput deletedInput =
        pipelineInputsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNull(deletedInput);
  }

  @Test
  // This test verifies that deleting pipeline inputs does NOT delete the associated pipeline run
  void deletingPipelineInputsDoesNotDeletePipelineRun() {
    // create and save a pipeline run
    UUID jobId = UUID.randomUUID();
    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(jobId);
    pipelineRun = pipelineRunsRepository.save(pipelineRun);

    // create and save inputs
    String inputsJson = "{\"input1\":\"gs://my-input-bucket/input1\",\"input2\":\"helloWorld\"}";
    PipelineInput input = new PipelineInput();
    input.setPipelineRunId(pipelineRun.getId());
    input.setInputs(inputsJson);
    pipelineInputsRepository.save(input);

    // verify pipeline run exists
    PipelineRun retrievedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNotNull(retrievedRun);

    // verify input exists
    PipelineInput retrievedInput =
        pipelineInputsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNotNull(retrievedInput);
    assertEquals(inputsJson, retrievedInput.getInputs());

    // delete input
    pipelineInputsRepository.delete(input);

    // flush and clear the persistence context
    entityManager.flush();
    entityManager.clear();

    // verify input is deleted
    PipelineInput deletedInput =
        pipelineInputsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNull(deletedInput);

    // verify pipeline run still exists
    PipelineRun nonDeletedPipelineRun =
        pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertNotNull(nonDeletedPipelineRun);
    assertEquals(jobId, nonDeletedPipelineRun.getJobId());
  }
}
