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

  private static final String GROOT_TEST_USER_ID = "i-am-groot";

  private QuotaAllocation newAllocation() {
    return new QuotaAllocation(
        PipelinesEnum.ARRAY_IMPUTATION,
        GROOT_TEST_USER_ID,
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
    assertEquals(GROOT_TEST_USER_ID, found.getUserId());
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
    QuotaAllocation defaultFreeAllocation = newAllocation();
    defaultFreeAllocation.setQuotaAllocated(2500);
    defaultFreeAllocation.setQuotaConsumed(1000);
    quotaAllocationsRepository.save(defaultFreeAllocation);

    QuotaAllocation paidAllocation =
        new QuotaAllocation(
            PipelinesEnum.ARRAY_IMPUTATION,
            GROOT_TEST_USER_ID,
            QuotaAllocationSourceEnum.EXTERNAL_PAID,
            500,
            200,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(paidAllocation);

    UserQuotaTotals totals =
        quotaAllocationsRepository.sumQuotaTotalsByUserIdAndPipelineName(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(3000, totals.totalAllocated());
    assertEquals(1200, totals.totalConsumed());
  }

  @Test
  void sumQuotaTotalsIsScopedToUserAndPipeline() {
    String otherUserId = "rocket-raccoon";

    QuotaAllocation firstUserArrayImputation = newAllocation();
    firstUserArrayImputation.setQuotaAllocated(2500);
    quotaAllocationsRepository.save(firstUserArrayImputation);

    // same user, different pipeline - should not be included
    QuotaAllocation firstUserLowPassImputation =
        new QuotaAllocation(
            PipelinesEnum.LOW_PASS_IMPUTATION,
            GROOT_TEST_USER_ID,
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
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(2500, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }

  @Test
  void findAllActiveAllocationsRequiresActiveTransactionForPessimisticLock() {
    // the query uses @Lock(PESSIMISTIC_WRITE), which requires an active transaction
    assertThrows(
        InvalidDataAccessApiUsageException.class,
        () ->
            quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
                GROOT_TEST_USER_ID,
                PipelinesEnum.ARRAY_IMPUTATION,
                QuotaAllocationStatusEnum.ACTIVE));
  }

  @Test
  @Transactional
  void findAllActiveAllocationsReturnsEmptyListOnNoMatch() {
    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
            "user-doesnt-exist", PipelinesEnum.ARRAY_IMPUTATION, QuotaAllocationStatusEnum.ACTIVE);

    assertTrue(result.isEmpty());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsReturnsOldestFirst() {
    Instant now = Instant.now();

    QuotaAllocation newest = newAllocation();
    newest.setQuotaSource(QuotaAllocationSourceEnum.EXTERNAL_PAID);
    quotaAllocationsRepository.save(newest);

    QuotaAllocation oldest = newAllocation();
    quotaAllocationsRepository.save(oldest);

    QuotaAllocation middle = newAllocation();
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
        quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION, QuotaAllocationStatusEnum.ACTIVE);

    // verify that the allocations are returned in chronological order
    assertEquals(3, result.size());
    assertEquals(oldest.getId(), result.get(0).getId());
    assertEquals(middle.getId(), result.get(1).getId());
    assertEquals(newest.getId(), result.get(2).getId());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsOrdersByCreatedWhileExcludingNonActive() {
    Instant now = Instant.now();

    QuotaAllocation newestActive = newAllocation();
    newestActive.setQuotaSource(QuotaAllocationSourceEnum.EXTERNAL_PAID);
    quotaAllocationsRepository.save(newestActive);

    // exhausted row sits chronologically between the two active rows
    QuotaAllocation exhausted = newAllocation();
    exhausted.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    exhausted.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    quotaAllocationsRepository.save(exhausted);

    QuotaAllocation oldestActive = newAllocation();
    quotaAllocationsRepository.save(oldestActive);

    setCreated(newestActive, now.plus(15, ChronoUnit.MINUTES));
    setCreated(exhausted, now.plus(10, ChronoUnit.MINUTES));
    setCreated(oldestActive, now);
    entityManager.flush();
    entityManager.clear();

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION, QuotaAllocationStatusEnum.ACTIVE);

    // exhausted row is excluded and the remaining rows are still in chronological order
    assertEquals(2, result.size());
    assertEquals(oldestActive.getId(), result.get(0).getId());
    assertEquals(newestActive.getId(), result.get(1).getId());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsFiltersByStatus() {
    QuotaAllocation activeAllocation = newAllocation();
    quotaAllocationsRepository.save(activeAllocation);

    QuotaAllocation exhaustedAllocation = newAllocation();
    exhaustedAllocation.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    exhaustedAllocation.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    quotaAllocationsRepository.save(exhaustedAllocation);

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION, QuotaAllocationStatusEnum.ACTIVE);

    assertEquals(1, result.size());
    assertEquals(activeAllocation.getId(), result.get(0).getId());
    assertEquals(QuotaAllocationStatusEnum.ACTIVE, result.get(0).getQuotaStatus());
  }

  @Test
  @Transactional
  void findAllActiveAllocationsFiltersByUserAndPipeline() {
    String otherUserId = "rocket-raccoon";

    QuotaAllocation grootArrayAllocation = newAllocation();
    quotaAllocationsRepository.save(grootArrayAllocation);

    // same user, different pipeline - should not be included
    QuotaAllocation grootLowPassAllocation =
        new QuotaAllocation(
            PipelinesEnum.LOW_PASS_IMPUTATION,
            GROOT_TEST_USER_ID,
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
        quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION, QuotaAllocationStatusEnum.ACTIVE);

    assertEquals(1, result.size());
    assertEquals(GROOT_TEST_USER_ID, result.get(0).getUserId());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, result.get(0).getPipelineName());
  }

  @Test
  void findAllAllocationsReturnsEmptyListOnNoMatch() {
    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameOrderByCreatedDesc(
            "not-groot-user", PipelinesEnum.ARRAY_IMPUTATION);

    assertTrue(result.isEmpty());
  }

  @Test
  @Transactional // required by setCreated helper method (native update)
  void findAllAllocationsReturnsNewestFirst() {
    Instant now = Instant.now();

    QuotaAllocation oldest = newAllocation();
    quotaAllocationsRepository.save(oldest);

    QuotaAllocation newest = newAllocation();
    newest.setQuotaSource(QuotaAllocationSourceEnum.EXTERNAL_PAID);
    quotaAllocationsRepository.save(newest);

    QuotaAllocation middle = newAllocation();
    middle.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    quotaAllocationsRepository.save(middle);

    // set created timestamps out of insertion order to verify the query sorts by created
    setCreated(oldest, now);
    setCreated(middle, now.plus(10, ChronoUnit.MINUTES));
    setCreated(newest, now.plus(15, ChronoUnit.MINUTES));
    entityManager.flush();
    entityManager.clear();

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameOrderByCreatedDesc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(3, result.size());
    assertEquals(newest.getId(), result.get(0).getId());
    assertEquals(middle.getId(), result.get(1).getId());
    assertEquals(oldest.getId(), result.get(2).getId());
  }

  @Test
  @Transactional // required by setCreated helper method (native update)
  void findAllAllocationsIncludesAllStatusesNewestFirst() {
    Instant now = Instant.now();

    QuotaAllocation active = newAllocation();
    quotaAllocationsRepository.save(active);

    QuotaAllocation exhausted = newAllocation();
    exhausted.setQuotaSource(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    exhausted.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    quotaAllocationsRepository.save(exhausted);

    setCreated(active, now);
    setCreated(exhausted, now.plus(5, ChronoUnit.MINUTES));
    entityManager.flush();
    entityManager.clear();

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameOrderByCreatedDesc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION);

    // non-active allocations are included
    assertEquals(2, result.size());
    assertEquals(exhausted.getId(), result.get(0).getId());
    assertEquals(QuotaAllocationStatusEnum.EXHAUSTED, result.get(0).getQuotaStatus());
    assertEquals(active.getId(), result.get(1).getId());
    assertEquals(QuotaAllocationStatusEnum.ACTIVE, result.get(1).getQuotaStatus());
  }

  @Test
  void findAllAllocationsFiltersByUserAndPipeline() {
    String otherUserId = "rocket-raccoon";

    QuotaAllocation grootArrayAllocation = newAllocation();
    quotaAllocationsRepository.save(grootArrayAllocation);

    // same user, different pipeline - should not be included
    QuotaAllocation grootLowPassAllocation =
        new QuotaAllocation(
            PipelinesEnum.LOW_PASS_IMPUTATION,
            GROOT_TEST_USER_ID,
            QuotaAllocationSourceEnum.DEFAULT_FREE,
            100,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(grootLowPassAllocation);

    // different user, same pipeline - should not be included
    QuotaAllocation raccoonArrayAllocation =
        new QuotaAllocation(
            PipelinesEnum.ARRAY_IMPUTATION,
            otherUserId,
            QuotaAllocationSourceEnum.DEFAULT_FREE,
            10000,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            null);
    quotaAllocationsRepository.save(raccoonArrayAllocation);

    List<QuotaAllocation> result =
        quotaAllocationsRepository.findByUserIdAndPipelineNameOrderByCreatedDesc(
            GROOT_TEST_USER_ID, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(1, result.size());
    assertEquals(grootArrayAllocation.getId(), result.get(0).getId());
    assertEquals(GROOT_TEST_USER_ID, result.get(0).getUserId());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, result.get(0).getPipelineName());
  }
}
