package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputsDefinition;
import java.util.List;
import org.springframework.data.repository.CrudRepository;

public interface PipelineInputsDefinitionsRepository
    extends CrudRepository<PipelineInputsDefinition, Long> {
  List<PipelineInputsDefinition> findAllByPipelineId(Pipeline pipeline);
}
