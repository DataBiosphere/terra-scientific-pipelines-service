package bio.terra.pipelines.db.repositories;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.QuotaAllocationDeduction;
import bio.terra.pipelines.db.entities.UserQuotaAllocation;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;

class QuotaAllocationDeductionsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private QuotaAllocationDeductionsRepository deductionsRepository;
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
                new QuotaAllocationDeduction(
                    allocationId, pipelineRun.getId(), 100, "run consumed quota"))
            .getId();

    QuotaAllocationDeduction found = deductionsRepository.findById(id).orElseThrow();
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
            .save(new QuotaAllocationDeduction(savedAllocationId(), null, 50, null))
            .getId();

    QuotaAllocationDeduction found = deductionsRepository.findById(id).orElseThrow();
    assertNull(found.getPipelineRunId());
    assertNull(found.getComments());
  }

  @Test
  void nonexistentAllocationIsRejected() {
    QuotaAllocationDeduction deduction = new QuotaAllocationDeduction(-1L, null, 50, null);
    assertThrows(DataIntegrityViolationException.class, () -> deductionsRepository.save(deduction));
  }

  @Test
  void nonexistentPipelineRunIsRejected() {
    QuotaAllocationDeduction deduction =
        new QuotaAllocationDeduction(savedAllocationId(), -1L, 50, null);
    assertThrows(DataIntegrityViolationException.class, () -> deductionsRepository.save(deduction));
  }
}
