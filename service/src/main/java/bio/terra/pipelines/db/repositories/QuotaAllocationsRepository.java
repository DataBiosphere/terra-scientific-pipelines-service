package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.model.UserQuotaTotals;
import jakarta.persistence.LockModeType;
import java.util.List;
import org.springframework.data.jpa.repository.Lock;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.CrudRepository;
import org.springframework.data.repository.query.Param;

public interface QuotaAllocationsRepository extends CrudRepository<QuotaAllocation, Long> {

  // CAST(... AS integer) is required: SUM() over an int-typed attribute resolves to Long at the
  // JPQL level, which would not match the QuotaTotals(int, int) constructor in the `new` expression
  // below without an explicit cast.
  @Query(
      "SELECT new bio.terra.pipelines.model.UserQuotaTotals("
          + "CAST(COALESCE(SUM(a.quotaAllocated), 0) AS integer), "
          + "CAST(COALESCE(SUM(a.quotaConsumed), 0) AS integer)) "
          + "FROM QuotaAllocation a WHERE a.userId = :userId AND a.pipelineName = :pipelineName")
  UserQuotaTotals sumQuotaTotalsByUserIdAndPipelineName(
      @Param("userId") String userId, @Param("pipelineName") PipelinesEnum pipelineName);

  // Locking is required to avoid race conditions when multiple threads are trying to
  // allocate quota for the same user and pipeline.
  // Note: Must be called within an active transaction for the lock to take effect
  @Lock(LockModeType.PESSIMISTIC_WRITE)
  List<QuotaAllocation> findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
      String userId, PipelinesEnum pipelineName, QuotaAllocationStatusEnum status);
}
