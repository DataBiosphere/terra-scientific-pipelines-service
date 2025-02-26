package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.CrudRepository;

public interface PipelineQuotasRepository extends CrudRepository<PipelineQuota, Long> {
  PipelineQuota findByPipelineName(PipelinesEnum pipelineName);

  @Query("SELECT p.quotaUnits FROM PipelineQuota p WHERE p.pipelineName = ?1")
  QuotaUnitsEnum findQuotaUnitsByPipeline(PipelinesEnum pipelineName);
}
