package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineInput;
import org.springframework.data.repository.CrudRepository;

public interface PipelineInputsRepository extends CrudRepository<PipelineInput, Long> {}
