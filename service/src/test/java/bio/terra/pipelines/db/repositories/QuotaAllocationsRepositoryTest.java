package bio.terra.pipelines.db.repositories;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.model.UserQuotaTotals;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import jakarta.persistence.EntityManager;
import java.time.Instant;
import java.time.temporal.ChronoUnit;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.InvalidDataAccessApiUsageException;
import org.springframework.transaction.annotation.Transactional;

class QuotaAllocationsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private QuotaAllocationsRepository quotaAllocationsRepository;

  @Autowired private EntityManager entityManager;

  private QuotaAllocation newAllocation() {
    return new QuotaAllocation(
        PipelinesEnum.ARRAY_IMPUTATION,
        "groot-user-id-123",
        QuotaAllocationSourceEnum.DEFAULT_FREE,
        2500,
        0,
        QuotaAllocationStatusEnum.ACTIVE,
        null);
  }

  // This method is used to set the created timestamp of a QuotaAllocation entity in the database.
  // This method allows for testing purposes to ensure that the created timestamps of different
  // allocations can be manipulated to test ordering and filtering logic. It uses a native
  // SQL query to update the created timestamp directly in the database, bypassing any
  // JPA entity management
  private void setCreated(QuotaAllocation allocation, Instant created) {
    entityManager
        .createNativeQuery("UPDATE quota_allocations SET created = :created WHERE id = :id")
        .setParameter("created", created)
        .setParameter("id", allocation.getId())
        .executeUpdate();
  }

  @Test
  void saveAllocation() {
    QuotaAllocation saved = quotaAllocationsRepository.save(newAllocation());
    assertNotNull(saved.getId());

    QuotaAllocation found = quotaAllocationsRepository.findById(saved.getId()).orElseThrow();
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, found.getPipelineName());
    assertEquals("groot-user-id-123", found.getUserId());
    assertEquals(QuotaAllocationSourceEnum.DEFAULT_FREE, found.getQuotaSource());
    assertEquals(2500, found.getQuotaAllocated());
    assertEquals(0, found.getQuotaConsumed());
    assertEquals(QuotaAllocationStatusEnum.ACTIVE, found.getQuotaStatus());
    assertNull(found.getComments());
    assertNotNull(found.getCreated());
    assertNotNull(found.getUpdated());
  }

  @Test
  void savingAllocationWithComments() {
    QuotaAllocation allocation = newAllocation();
    allocation.setComments("I am Groot");

    Long id = quotaAllocationsRepository.save(allocation).getId();

    assertEquals("I am Groot", quotaAllocationsRepository.findById(id).orElseThrow().getComments());
  }

  @Test
  void allQuotaSourcesCanBePersisted() {
    for (QuotaAllocationSourceEnum source : QuotaAllocationSourceEnum.values()) {
      QuotaAllocation allocation = newAllocation();
      allocation.setQuotaSource(source);
      Long id = quotaAllocationsRepository.save(allocation).getId();
      assertEquals(source, quotaAllocationsRepository.findById(id).orElseThrow().getQuotaSource());
    }
  }

  @Test
  void updatingRowUpdatesUpdatedTimestamp() {
    Long id = quotaAllocationsRepository.save(newAllocation()).getId();
    QuotaAllocation quotaAllocation = quotaAllocationsRepository.findById(id).orElseThrow();
    Instant createdTimestamp = quotaAllocation.getCreated();
    Instant updatedTimestampBeforeUpdate = quotaAllocation.getUpdated();

    quotaAllocation.setQuotaConsumed(2500);
    quotaAllocation.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    quotaAllocationsRepository.save(quotaAllocation);

    QuotaAllocation quotaAllocationAfterUpdated =
        quotaAllocationsRepository.findById(id).orElseThrow();
    assertEquals(2500, quotaAllocationAfterUpdated.getQuotaConsumed());
    assertEquals(QuotaAllocationStatusEnum.EXHAUSTED, quotaAllocationAfterUpdated.getQuotaStatus());
    assertTrue(quotaAllocationAfterUpdated.getUpdated().isAfter(updatedTimestampBeforeUpdate));
    assertEquals(createdTimestamp, quotaAllocationAfterUpdated.getCreated());
  }

  @Test
  void sumQuotaTotalsReturnsZeroWhenNoAllocationsExist() {
    UserQuotaTotals totals =
        quotaAllocationsRepository.sumQuotaTotalsByUserIdAndPipelineName(
            "groot-has-no-quota", PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(0, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }

  @Test
  void sumQuotaTotalsAggregatesAcrossMultipleAllocations() {
    String userId = "i-am-groot";

    QuotaAllocation defaultFreeAllocation = newAllocation();
    defaultFreeAllocation.setUserId(userId);
    defaultFreeAllocation.setQuotaAllocated(2500);
    defaultFreeAllocation.setQuotaConsumed(1000);
    quotaAllocationsRepository.save(defaultFreeAllocation);

    QuotaAllocation paidAllocation =
        new QuotaAllocation(
            PipelinesEnum.ARRAY_IMPUTATION,
            userId,
            QuotaAllocationSourceEnum.EXTERNAL_PAID,
            500,
            200,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(paidAllocation);

    UserQuotaTotals totals =
        quotaAllocationsRepository.sumQuotaTotalsByUserIdAndPipelineName(
            userId, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(3000, totals.totalAllocated());
    assertEquals(1200, totals.totalConsumed());
  }

  @Test
  void sumQuotaTotalsIsScopedToUserAndPipeline() {
    String userId = "i-am-groot";
    String otherUserId = "rocket-raccoon";

    QuotaAllocation firstUserArrayImputation = newAllocation();
    firstUserArrayImputation.setUserId(userId);
    firstUserArrayImputation.setQuotaAllocated(2500);
    quotaAllocationsRepository.save(firstUserArrayImputation);

    // same user, different pipeline - should not be included
    QuotaAllocation firstUserLowPassImputation =
        new QuotaAllocation(
            PipelinesEnum.LOW_PASS_IMPUTATION,
            userId,
            QuotaAllocationSourceEnum.DEFAULT_FREE,
            100,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(firstUserLowPassImputation);

    // different user, same pipeline - should not be included
    QuotaAllocation forOtherUser =
        new QuotaAllocation(
            PipelinesEnum.ARRAY_IMPUTATION,
            otherUserId,
            QuotaAllocationSourceEnum.DEFAULT_FREE,
            10000,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(forOtherUser);

    UserQuotaTotals totals =
        quotaAllocationsRepository.sumQuotaTotalsByUserIdAndPipelineName(
            userId, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(2500, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }

  @Test
  void findAllActiveAllocationsRequiresActiveTransactionForPessimisticLock() {
    // the query uses @Lock(PESSIMISTIC_WRITE), which requires an active transaction
    assertThrows(
        InvalidDataAccessApiUsageException.class,
        () ->
            quotaAllocationsRepository
                .findActiveAllocationsByUserIdAndPipelineNameOrderByCreatedAsc(
                    "i-am-groot", PipelinesEnum.ARRAY_IMPUTATION));
  }

  @Test
  @Transactional
  void findAllActiveAllocationsReturnsEmptyListOnNoMatch() {
    List<QuotaAllocation> result =
        quotaAllocationsRepository.findActiveAllocationsByUserIdAndPipelineNameOrderByCreatedAsc(
            "user-doesnt-exist", PipelinesEnum.ARRAY_IMPUTATION);

    assertTrue(result.isEmpty());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsReturnsOldestFirst() {
    String userId = "i-am-groot";
    Instant now = Instant.now();

    QuotaAllocation newest = newAllocation();
    newest.setUserId(userId);
    newest.setQuotaSource(QuotaAllocationSourceEnum.EXTERNAL_PAID);
    quotaAllocationsRepository.save(newest);

    QuotaAllocation oldest = newAllocation();
    oldest.setUserId(userId);
    quotaAllocationsRepository.save(oldest);

    QuotaAllocation middle = newAllocation();
    middle.setUserId(userId);
    middle.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    quotaAllocationsRepository.save(middle);

    // update the created timestamps to ensure they are out of order. This is to test
    // that the query returns them in chronological order regardless of the order they were saved in
    setCreated(newest, now.plus(15, ChronoUnit.MINUTES));
    setCreated(oldest, now);
    setCreated(middle, now.plus(10, ChronoUnit.MINUTES));
    entityManager.flush();
    entityManager.clear();

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findActiveAllocationsByUserIdAndPipelineNameOrderByCreatedAsc(
            userId, PipelinesEnum.ARRAY_IMPUTATION);

    // verify that the allocations are returned in chronological order
    assertEquals(3, result.size());
    assertEquals(oldest.getId(), result.get(0).getId());
    assertEquals(middle.getId(), result.get(1).getId());
    assertEquals(newest.getId(), result.get(2).getId());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsOrdersByCreatedWhileExcludingNonActive() {
    String userId = "i-am-groot";
    Instant now = Instant.now();

    QuotaAllocation newestActive = newAllocation();
    newestActive.setUserId(userId);
    newestActive.setQuotaSource(QuotaAllocationSourceEnum.EXTERNAL_PAID);
    quotaAllocationsRepository.save(newestActive);

    // exhausted row sits chronologically between the two active rows
    QuotaAllocation exhausted = newAllocation();
    exhausted.setUserId(userId);
    exhausted.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    exhausted.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    quotaAllocationsRepository.save(exhausted);

    QuotaAllocation oldestActive = newAllocation();
    oldestActive.setUserId(userId);
    quotaAllocationsRepository.save(oldestActive);

    setCreated(newestActive, now.plus(15, ChronoUnit.MINUTES));
    setCreated(exhausted, now.plus(10, ChronoUnit.MINUTES));
    setCreated(oldestActive, now);
    entityManager.flush();
    entityManager.clear();

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findActiveAllocationsByUserIdAndPipelineNameOrderByCreatedAsc(
            userId, PipelinesEnum.ARRAY_IMPUTATION);

    // exhausted row is excluded and the remaining rows are still in chronological order
    assertEquals(2, result.size());
    assertEquals(oldestActive.getId(), result.get(0).getId());
    assertEquals(newestActive.getId(), result.get(1).getId());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsFiltersByStatus() {
    String userId = "i-am-groot";

    QuotaAllocation activeAllocation = newAllocation();
    activeAllocation.setUserId(userId);
    quotaAllocationsRepository.save(activeAllocation);

    QuotaAllocation exhaustedAllocation = newAllocation();
    exhaustedAllocation.setUserId(userId);
    exhaustedAllocation.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    exhaustedAllocation.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    quotaAllocationsRepository.save(exhaustedAllocation);

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findActiveAllocationsByUserIdAndPipelineNameOrderByCreatedAsc(
            userId, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(1, result.size());
    assertEquals(activeAllocation.getId(), result.get(0).getId());
    assertEquals(QuotaAllocationStatusEnum.ACTIVE, result.get(0).getQuotaStatus());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsFiltersByUserAndPipeline() {
    String userId = "i-am-groot";
    String otherUserId = "rocket-raccoon";

    QuotaAllocation grootArrayAllocation = newAllocation();
    grootArrayAllocation.setUserId(userId);
    quotaAllocationsRepository.save(grootArrayAllocation);

    // same user, different pipeline - should not be included
    QuotaAllocation grootLowPassAllocation =
        new QuotaAllocation(
            PipelinesEnum.LOW_PASS_IMPUTATION,
            userId,
            QuotaAllocationSourceEnum.DEFAULT_FREE,
            100,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(grootLowPassAllocation);

    // different user, same pipeline - should not be included
    QuotaAllocation racoonArrayAllocation =
        new QuotaAllocation(
            PipelinesEnum.ARRAY_IMPUTATION,
            otherUserId,
            QuotaAllocationSourceEnum.DEFAULT_FREE,
            10000,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(racoonArrayAllocation);

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findActiveAllocationsByUserIdAndPipelineNameOrderByCreatedAsc(
            userId, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(1, result.size());
    assertEquals(userId, result.get(0).getUserId());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, result.get(0).getPipelineName());
  }
}
