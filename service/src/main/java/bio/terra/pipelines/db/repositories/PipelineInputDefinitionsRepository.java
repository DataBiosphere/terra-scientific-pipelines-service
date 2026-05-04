package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import org.springframework.data.repository.CrudRepository;
import org.springframework.transaction.annotation.Transactional;

public interface PipelineInputDefinitionsRepository
    extends CrudRepository<PipelineInputDefinition, Long> {
  @Transactional
  void deleteByPipelineId(Long pipelineId);
}
