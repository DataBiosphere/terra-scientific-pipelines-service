package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.model.UserQuotaTotals;
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
}
