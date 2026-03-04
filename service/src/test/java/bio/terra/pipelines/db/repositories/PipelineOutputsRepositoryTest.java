package bio.terra.pipelines.db.repositories;

import static org.assertj.core.api.Assertions.assertThat;

import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import jakarta.persistence.EntityManager;
import java.util.List;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.transaction.annotation.Transactional;

@Transactional
class PipelineOutputsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Autowired private PipelineOutputsRepository pipelineOutputsRepository;
  @Autowired private EntityManager entityManager;

  @Test
  // This test verifies that when a pipeline run is deleted, all associated outputs are also deleted
  // due to the cascade delete configuration
  void deletingPipelineRunDeletesAssociatedOutputs() {
    // create and save a pipeline pipelineRun
    UUID jobId = UUID.randomUUID();
    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(jobId);
    pipelineRun = pipelineRunsRepository.save(pipelineRun);

    // create and save outputs
    PipelineOutput output1 = new PipelineOutput();
    output1.setPipelineRunId(pipelineRun.getId());
    output1.setOutputName("output1");
    output1.setOutputValue("gs://my-output-bucket/output1");
    PipelineOutput output2 = new PipelineOutput();
    output2.setPipelineRunId(pipelineRun.getId());
    output2.setOutputName("output2");
    output2.setOutputValue("helloWorld");
    pipelineOutputsRepository.saveAll(List.of(output1, output2));

    // verify pipeline pipelineRun exists
    PipelineRun retrievedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(retrievedRun).isNotNull();
    assertThat(retrievedRun.getJobId()).isEqualTo(jobId);

    // verify outputs exist
    List<PipelineOutput> outputs =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(pipelineRun.getId());
    assertThat(outputs).hasSize(2);
    assertThat(outputs)
        .extracting(PipelineOutput::getOutputName)
        .containsExactlyInAnyOrder("output1", "output2");
    assertThat(outputs)
        .extracting(PipelineOutput::getOutputValue)
        .containsExactlyInAnyOrder("gs://my-output-bucket/output1", "helloWorld");

    // delete pipeline pipelineRun
    pipelineRunsRepository.delete(pipelineRun);

    // flush and clear the persistence context to ensure we are fetching fresh data from the
    // database
    entityManager.flush();
    entityManager.clear();

    // verify pipeline pipelineRun is deleted
    PipelineRun deletedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(deletedRun).isNull();

    // verify cascade delete works and outputs are deleted
    List<PipelineOutput> deletedOutputs =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(pipelineRun.getId());
    assertThat(deletedOutputs).isEmpty();
  }

  @Test
  // This test verifies that deleting pipeline outputs does NOT delete the associated pipeline run
  void deletingPipelineOutputsDoesNotDeletePipelineRun() {
    // create and save a pipeline run
    UUID jobId = UUID.randomUUID();
    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(jobId);
    pipelineRun = pipelineRunsRepository.save(pipelineRun);

    // create and save outputs
    PipelineOutput output1 = new PipelineOutput();
    output1.setPipelineRunId(pipelineRun.getId());
    output1.setOutputName("output1");
    output1.setOutputValue("gs://my-output-bucket/output1");
    PipelineOutput output2 = new PipelineOutput();
    output2.setPipelineRunId(pipelineRun.getId());
    output2.setOutputName("output2");
    output2.setOutputValue("helloWorld");
    pipelineOutputsRepository.saveAll(List.of(output1, output2));

    // verify pipeline pipelineRun exists
    PipelineRun retrievedRun = pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(retrievedRun).isNotNull();
    assertThat(retrievedRun.getJobId()).isEqualTo(jobId);

    // verify outputs exist
    List<PipelineOutput> outputs =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(pipelineRun.getId());
    assertThat(outputs).hasSize(2);

    // delete all outputs
    pipelineOutputsRepository.deleteAll(outputs);

    // flush and clear the persistence context to ensure we are fetching fresh data from the
    // database
    entityManager.flush();
    entityManager.clear();

    // verify outputs are deleted
    List<PipelineOutput> deletedOutputs =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(pipelineRun.getId());
    assertThat(deletedOutputs).isEmpty();

    // verify pipeline run still exists
    PipelineRun nonDeletedPipelineRun =
        pipelineRunsRepository.findById(pipelineRun.getId()).orElse(null);
    assertThat(nonDeletedPipelineRun).isNotNull();
    assertThat(nonDeletedPipelineRun.getJobId()).isEqualTo(jobId);
  }
}
