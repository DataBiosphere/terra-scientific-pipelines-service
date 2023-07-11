package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineInput;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface PipelineInputsRepository extends CrudRepository<PipelineInput, UUID> {}
