package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import org.springframework.data.repository.CrudRepository;

public interface PipelineInputDefinitionsRepository
    extends CrudRepository<PipelineInputDefinition, Long> {}
