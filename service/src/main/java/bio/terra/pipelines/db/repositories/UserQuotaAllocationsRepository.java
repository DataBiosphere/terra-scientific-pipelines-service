package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.db.entities.UserQuotaAllocation;
import bio.terra.pipelines.model.UserQuotaTotals;
import java.util.Optional;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.CrudRepository;
import org.springframework.data.repository.query.Param;

public interface UserQuotaAllocationsRepository extends CrudRepository<UserQuotaAllocation, Long> {

  Optional<UserQuotaAllocation> findByUserIdAndPipelineNameAndQuotaSource(
      String userId, PipelinesEnum pipelineName, QuotaAllocationSourceEnum quotaSource);

  // CAST(... AS integer) is required: SUM() over an int-typed attribute resolves to Long at the
  // JPQL level, which would not match the QuotaTotals(int, int) constructor in the `new` expression
  // below without an explicit cast.
  @Query(
      "SELECT new bio.terra.pipelines.model.UserQuotaTotals("
          + "CAST(COALESCE(SUM(a.quotaAllocated), 0) AS integer), "
          + "CAST(COALESCE(SUM(a.quotaConsumed), 0) AS integer)) "
          + "FROM UserQuotaAllocation a WHERE a.userId = :userId AND a.pipelineName = :pipelineName")
  UserQuotaTotals sumQuotaTotalsByUserIdAndPipelineName(
      @Param("userId") String userId, @Param("pipelineName") PipelinesEnum pipelineName);
}
