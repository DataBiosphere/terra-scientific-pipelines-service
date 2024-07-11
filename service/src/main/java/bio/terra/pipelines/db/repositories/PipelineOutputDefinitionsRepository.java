package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import org.springframework.data.repository.CrudRepository;

public interface PipelineOutputDefinitionsRepository
    extends CrudRepository<PipelineOutputDefinition, Long> {}
