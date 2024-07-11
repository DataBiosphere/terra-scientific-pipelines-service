package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineInputs;
import org.springframework.data.repository.CrudRepository;

public interface PipelineInputsRepository extends CrudRepository<PipelineInputs, Long> {}
