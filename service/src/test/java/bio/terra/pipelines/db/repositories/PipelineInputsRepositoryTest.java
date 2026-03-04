package bio.terra.pipelines.db.repositories;

import static org.assertj.core.api.Assertions.assertThat;

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
public class PipelineInputsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Autowired private PipelineInputsRepository pipelineInputsRepository;
  @Autowired private EntityManager entityManager;

  @Test
  // This test verifies that when a pipeline run is deleted, all associated inputs are also deleted
  // due to the cascade delete configuration
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
    assertThat(retrievedRun).isNotNull();

    // verify input exists
    PipelineInput retrievedInput =
        pipelineInputsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(retrievedInput).isNotNull();
    assertThat(retrievedRun.getJobId()).isEqualTo(jobId);
    assertThat(retrievedInput.getInputs()).isEqualTo(inputsJson);

    // delete pipeline run
    pipelineRunsRepository.delete(pipelineRun);

    // flush and clear the persistence context to ensure we are fetching fresh data from the
    // database
    entityManager.flush();
    entityManager.clear();

    // verify pipeline run is deleted
    PipelineRun deletedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(deletedRun).isNull();

    // verify cascade delete - input should also be deleted
    PipelineInput deletedInput =
        pipelineInputsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(deletedInput).isNull();
  }
}
