package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.QuotaAllocationDeduction;
import org.springframework.data.repository.CrudRepository;

public interface QuotaAllocationDeductionsRepository
    extends CrudRepository<QuotaAllocationDeduction, Long> {}
