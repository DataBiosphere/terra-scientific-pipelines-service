package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.QuotaAllocationDeduction;
import java.util.List;
import org.springframework.data.repository.CrudRepository;

public interface QuotaAllocationDeductionsRepository
    extends CrudRepository<QuotaAllocationDeduction, Long> {

  boolean existsByPipelineRunId(Long pipelineRunId);

  List<QuotaAllocationDeduction> findByPipelineRunId(Long pipelineRunId);
}
