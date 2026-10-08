package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.UserQuotaAllocationDeduction;
import org.springframework.data.repository.CrudRepository;

public interface UserQuotaAllocationDeductionsRepository
    extends CrudRepository<UserQuotaAllocationDeduction, Long> {}
