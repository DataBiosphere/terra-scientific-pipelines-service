package bio.terra.pipelines.db.repositories;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.db.entities.UserQuotaAllocation;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.time.Instant;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class UserQuotaAllocationsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private UserQuotaAllocationsRepository userQuotaAllocationsRepository;

  private UserQuotaAllocation newAllocation() {
    return new UserQuotaAllocation(
        PipelinesEnum.ARRAY_IMPUTATION,
        "groot-user-id-123",
        QuotaAllocationSourceEnum.DEFAULT_FREE,
        2500,
        0,
        QuotaAllocationStatusEnum.ACTIVE,
        null);
  }

  @Test
  void saveAllocation() {
    UserQuotaAllocation saved = userQuotaAllocationsRepository.save(newAllocation());
    assertNotNull(saved.getId());

    UserQuotaAllocation found =
        userQuotaAllocationsRepository.findById(saved.getId()).orElseThrow();
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
    UserQuotaAllocation allocation = newAllocation();
    allocation.setComments("I am Groot");

    Long id = userQuotaAllocationsRepository.save(allocation).getId();

    assertEquals(
        "I am Groot", userQuotaAllocationsRepository.findById(id).orElseThrow().getComments());
  }

  @Test
  void allQuotaSourcesCanBePersisted() {
    for (QuotaAllocationSourceEnum source : QuotaAllocationSourceEnum.values()) {
      UserQuotaAllocation allocation = newAllocation();
      allocation.setQuotaSource(source);
      Long id = userQuotaAllocationsRepository.save(allocation).getId();
      assertEquals(
          source, userQuotaAllocationsRepository.findById(id).orElseThrow().getQuotaSource());
    }
  }

  @Test
  void updatingRowUpdatesUpdatedTimestamp() {
    Long id = userQuotaAllocationsRepository.save(newAllocation()).getId();
    UserQuotaAllocation quotaAllocation = userQuotaAllocationsRepository.findById(id).orElseThrow();
    Instant createdTimestamp = quotaAllocation.getCreated();
    Instant updatedTimestampBeforeUpdate = quotaAllocation.getUpdated();

    quotaAllocation.setQuotaConsumed(2500);
    quotaAllocation.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    userQuotaAllocationsRepository.save(quotaAllocation);

    UserQuotaAllocation quotaAllocationAfterUpdated =
        userQuotaAllocationsRepository.findById(id).orElseThrow();
    assertEquals(2500, quotaAllocationAfterUpdated.getQuotaConsumed());
    assertEquals(QuotaAllocationStatusEnum.EXHAUSTED, quotaAllocationAfterUpdated.getQuotaStatus());
    assertTrue(quotaAllocationAfterUpdated.getUpdated().isAfter(updatedTimestampBeforeUpdate));
    assertEquals(createdTimestamp, quotaAllocationAfterUpdated.getCreated());
  }
}
