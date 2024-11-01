package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import org.springframework.data.repository.CrudRepository;

public interface PipelineQuotasRepository extends CrudRepository<PipelineQuota, Long> {
  PipelineQuota findByPipelineName(PipelinesEnum pipelineName);
}
