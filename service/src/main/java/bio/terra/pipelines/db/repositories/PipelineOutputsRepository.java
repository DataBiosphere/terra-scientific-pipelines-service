package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineOutput;
import org.springframework.data.repository.CrudRepository;

public interface PipelineOutputsRepository extends CrudRepository<PipelineOutput, Long> {
  PipelineOutput findPipelineOutputsByJobId(Long jobId);
}
