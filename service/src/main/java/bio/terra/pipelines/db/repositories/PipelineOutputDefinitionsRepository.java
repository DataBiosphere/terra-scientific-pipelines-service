package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import org.springframework.data.repository.CrudRepository;
import org.springframework.transaction.annotation.Transactional;

public interface PipelineOutputDefinitionsRepository
    extends CrudRepository<PipelineOutputDefinition, Long> {
  @Transactional
  void deleteByPipelineId(Long pipelineId);
}
