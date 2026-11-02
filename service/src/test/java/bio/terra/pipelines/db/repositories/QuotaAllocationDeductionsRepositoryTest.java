package bio.terra.pipelines.db.repositories;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.db.entities.QuotaAllocationDeduction;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;

class QuotaAllocationDeductionsRepositoryTest extends BaseEmbeddedDbTest {

  @Autowired private QuotaAllocationDeductionsRepository deductionsRepository;
  @Autowired private QuotaAllocationsRepository allocationsRepository;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;

  private Long savedAllocationId(QuotaAllocationSourceEnum source) {
    return allocationsRepository
        .save(
            new QuotaAllocation(
                PipelinesEnum.ARRAY_IMPUTATION,
                "groot-user-id-123",
                source,
                2500,
                0,
                QuotaAllocationStatusEnum.ACTIVE,
                null))
        .getId();
  }

  @Test
  void saveDeductionWithPipelineRun() {
    Long allocationId = savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE);
    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));

    Long id =
        deductionsRepository
            .save(
                new QuotaAllocationDeduction(
                    allocationId, pipelineRun.getId(), 100, "run consumed quota"))
            .getId();

    QuotaAllocationDeduction found = deductionsRepository.findById(id).orElseThrow();
    assertEquals(allocationId, found.getQuotaAllocationId());
    assertEquals(pipelineRun.getId(), found.getPipelineRunId());
    assertEquals(100, found.getAmount());
    assertEquals("run consumed quota", found.getComments());
    assertNotNull(found.getCreated());
  }

  @Test
  void saveDeductionWithoutPipelineRun() {
    Long id =
        deductionsRepository
            .save(
                new QuotaAllocationDeduction(
                    savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE), null, 50, null))
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
        new QuotaAllocationDeduction(
            savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE), -1L, 50, null);
    assertThrows(DataIntegrityViolationException.class, () -> deductionsRepository.save(deduction));
  }

  @Test
  void existsByPipelineRunIdReturnsTrueWhenDeductionsExist() {
    Long allocationId = savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE);
    Long secondAllocationId = savedAllocationId(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);
    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));
    deductionsRepository.save(
        new QuotaAllocationDeduction(allocationId, pipelineRun.getId(), 100, null));
    deductionsRepository.save(
        new QuotaAllocationDeduction(secondAllocationId, pipelineRun.getId(), 50, null));

    assertTrue(deductionsRepository.existsByPipelineRunId(pipelineRun.getId()));
  }

  @Test
  void existsByPipelineRunIdReturnsFalseWhenPipelineRunHasNoDeductions() {
    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));

    assertFalse(deductionsRepository.existsByPipelineRunId(pipelineRun.getId()));
  }

  @Test
  void existsByPipelineRunIdReturnsFalseWhenOnlyOtherRunsHaveDeductions() {
    Long allocationId = savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE);
    PipelineRun runWithDeduction =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));
    PipelineRun runWithoutDeduction =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));
    deductionsRepository.save(
        new QuotaAllocationDeduction(allocationId, runWithDeduction.getId(), 100, null));

    assertTrue(deductionsRepository.existsByPipelineRunId(runWithDeduction.getId()));
    assertFalse(deductionsRepository.existsByPipelineRunId(runWithoutDeduction.getId()));
  }

  @Test
  void existsByPipelineRunIdIgnoresDeductionsWithoutPipelineRun() {
    // e.g. deductions recorded by an admin are not tied to any pipeline run
    deductionsRepository.save(
        new QuotaAllocationDeduction(
            savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE),
            null,
            50,
            "manual admin deduction"));
    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));

    assertFalse(deductionsRepository.existsByPipelineRunId(pipelineRun.getId()));
  }

  @Test
  void findByPipelineRunIdReturnsEmptyListWhenRunHasNoDeductions() {
    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));

    assertTrue(deductionsRepository.findByPipelineRunId(pipelineRun.getId()).isEmpty());
  }

  @Test
  void findByPipelineRunIdReturnsAllDeductionsForRun() {
    Long allocationId = savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE);
    Long secondAllocationId = savedAllocationId(QuotaAllocationSourceEnum.DEVELOPMENT_FREE);

    PipelineRun pipelineRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));
    deductionsRepository.save(
        new QuotaAllocationDeduction(allocationId, pipelineRun.getId(), 100, null));
    deductionsRepository.save(
        new QuotaAllocationDeduction(secondAllocationId, pipelineRun.getId(), 50, null));

    List<QuotaAllocationDeduction> result =
        deductionsRepository.findByPipelineRunId(pipelineRun.getId());

    assertEquals(2, result.size());
    assertTrue(result.stream().allMatch(d -> pipelineRun.getId().equals(d.getPipelineRunId())));
    assertEquals(150, result.stream().mapToInt(QuotaAllocationDeduction::getAmount).sum());
  }

  @Test
  void findByPipelineRunIdOnlyReturnsDeductionsForThatRun() {
    Long allocationId = savedAllocationId(QuotaAllocationSourceEnum.DEFAULT_FREE);

    PipelineRun firstRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));
    PipelineRun secondRun =
        pipelineRunsRepository.save(TestUtils.createNewPipelineRunWithJobId(UUID.randomUUID()));

    deductionsRepository.save(
        new QuotaAllocationDeduction(allocationId, firstRun.getId(), 100, null));
    deductionsRepository.save(
        new QuotaAllocationDeduction(allocationId, secondRun.getId(), 200, null));

    // a deduction that isn't tied to any run, e.g. recorded by an admin
    deductionsRepository.save(new QuotaAllocationDeduction(allocationId, null, 300, "admin"));

    List<QuotaAllocationDeduction> firstRunDeductions =
        deductionsRepository.findByPipelineRunId(firstRun.getId());
    assertEquals(1, firstRunDeductions.size());
    assertEquals(100, firstRunDeductions.get(0).getAmount());

    List<QuotaAllocationDeduction> secondRunDeductions =
        deductionsRepository.findByPipelineRunId(secondRun.getId());
    assertEquals(1, secondRunDeductions.size());
    assertEquals(200, secondRunDeductions.get(0).getAmount());
  }
}
