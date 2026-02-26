package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineOutput;
import java.util.List;
import org.springframework.data.repository.CrudRepository;

public interface PipelineOutputsRepository extends CrudRepository<PipelineOutput, Long> {

  List<PipelineOutput> findPipelineOutputsByPipelineRunsId(Long pipelineRunsId);
}
