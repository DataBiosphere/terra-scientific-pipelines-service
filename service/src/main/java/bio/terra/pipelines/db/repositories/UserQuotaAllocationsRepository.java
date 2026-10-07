package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.UserQuotaAllocation;
import org.springframework.data.repository.CrudRepository;

public interface UserQuotaAllocationsRepository extends CrudRepository<UserQuotaAllocation, Long> {}
