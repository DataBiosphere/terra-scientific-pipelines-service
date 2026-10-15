package bio.terra.pipelines.db.repositories;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.UserQuotaAllocation;
import bio.terra.pipelines.db.entities.UserQuotaAllocationDeduction;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;

class UserQuotaAllocationDeductionsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private UserQuotaAllocationDeductionsRepository deductionsRepository;
  @Autowired private UserQuotaAllocationsRepository allocationsRepository;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;

  private Long savedAllocationId() {
    return allocationsRepository
        .save(
            new UserQuotaAllocation(
                PipelinesEnum.ARRAY_IMPUTATION,
                "groot-user-id-123",
                QuotaAllocationSourceEnum.DEFAULT_FREE,
                2500,
                0,
                QuotaAllocationStatusEnum.ACTIVE,
                null))
        .getId();
  }

  @Test
  void saveDeductionWithPipelineRun() {
    Long allocationId = savedAllocationId();
    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));

    Long id =
        deductionsRepository
            .save(
                new UserQuotaAllocationDeduction(
                    allocationId, pipelineRun.getId(), 100, "run consumed quota"))
            .getId();

    UserQuotaAllocationDeduction found = deductionsRepository.findById(id).orElseThrow();
    assertEquals(allocationId, found.getUserQuotaAllocationId());
    assertEquals(pipelineRun.getId(), found.getPipelineRunId());
    assertEquals(100, found.getAmount());
    assertEquals("run consumed quota", found.getComments());
    assertNotNull(found.getCreated());
  }

  @Test
  void saveDeductionWithoutPipelineRun() {
    Long id =
        deductionsRepository
            .save(new UserQuotaAllocationDeduction(savedAllocationId(), null, 50, null))
            .getId();

    UserQuotaAllocationDeduction found = deductionsRepository.findById(id).orElseThrow();
    assertNull(found.getPipelineRunId());
    assertNull(found.getComments());
  }

  @Test
  void nonexistentAllocationIsRejected() {
    UserQuotaAllocationDeduction deduction = new UserQuotaAllocationDeduction(-1L, null, 50, null);
    assertThrows(DataIntegrityViolationException.class, () -> deductionsRepository.save(deduction));
  }

  @Test
  void nonexistentPipelineRunIsRejected() {
    UserQuotaAllocationDeduction deduction =
        new UserQuotaAllocationDeduction(savedAllocationId(), -1L, 50, null);
    assertThrows(DataIntegrityViolationException.class, () -> deductionsRepository.save(deduction));
  }
}
