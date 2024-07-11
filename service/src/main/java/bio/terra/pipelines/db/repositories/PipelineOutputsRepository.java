package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineOutputs;
import org.springframework.data.repository.CrudRepository;

public interface PipelineOutputsRepository extends CrudRepository<PipelineOutputs, Long> {
  PipelineOutputs findPipelineOutputsByJobId(Long jobId);
}
