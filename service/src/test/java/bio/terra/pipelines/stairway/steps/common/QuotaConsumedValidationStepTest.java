package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.testutils.TestUtils.createNewPipelineRunWithJobId;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.mockito.Mockito.any;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.model.PipelineQuota;
import bio.terra.pipelines.model.QuotaConsumptionRecord;
import bio.terra.pipelines.model.UserQuotaTotals;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.*;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.CannotAcquireLockException;

class QuotaConsumedValidationStepTest extends BaseEmbeddedDbTest {

  @Autowired private QuotasService quotasService;
  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;

  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private void putQuotaOutputs(String quotaConsumed) {
    flightContext
        .getWorkingMap()
        .put(
            WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS,
            new HashMap<>(Map.of("quotaConsumed", quotaConsumed)));
  }

  // a QuotasService mock where the user has 2500 quota with none consumed, so that the step passes
  // its quota checks and the test can control what happens when the quota is consumed
  private QuotasService mockQuotasServiceWithQuotaAvailable() {
    QuotasService mockQuotasService = mock(QuotasService.class);
    when(mockQuotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION))
        .thenReturn(quotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION));
    when(mockQuotasService.getQuotaTotalsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION))
        .thenReturn(new UserQuotaTotals(2500, 0));
    return mockQuotasService;
  }

  private void assertNothingStoredInWorkingMap() {
    FlightMap workingMap = flightContext.getWorkingMap();
    assertFalse(workingMap.containsKey(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS));
    assertFalse(workingMap.containsKey(WdlBasedPipelineJobMapKeys.RAW_QUOTA_CONSUMED));
  }

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    // set up the pipeline run
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    // the step doesn't create a default allocation, so make sure the user has one
    quotasService.getOrInitializeQuotaTotalsForUserAndPipeline(
        TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
  }

  @Test
  void doStepSuccessUsePipelineQuotaFloor() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quotaConsumed to 30 which is below the pipeline min quota of 500
    final Integer quotaConsumedBelowPipelineMinQuota = 30;
    final Map<String, String> quotaOutputs =
        new HashMap<>(Map.of("quotaConsumed", quotaConsumedBelowPipelineMinQuota.toString()));
    flightContext.getWorkingMap().put(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // before running make sure quota consumed for user is 0
    UserQuotaTotals userQuotaTotals =
        quotasService.getQuotaTotalsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(0, userQuotaTotals.totalConsumed());

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // after running make sure quota for user is 500 from the pipeline min quota
    userQuotaTotals =
        quotasService.getQuotaTotalsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
    PipelineQuota pipelineQuota = quotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(pipelineQuota.getMinQuotaConsumed(), userQuotaTotals.totalConsumed());

    // make sure the raw quota assumed was saved for the pipeline
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_1_ID)
            .orElseThrow();
    assertEquals(quotaConsumedBelowPipelineMinQuota, updatedPipelineRun.getRawQuotaConsumed());

    assertEquals(
        quotaConsumedBelowPipelineMinQuota,
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.RAW_QUOTA_CONSUMED, Integer.class));

    // assert that the quota consumption records are stored in the working map
    List<QuotaConsumptionRecord> records =
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, List.class);

    assertEquals(1, records.size());

    QuotaConsumptionRecord record = records.get(0);
    assertEquals(500, record.amountConsumed());
  }

  @Test
  void doStepSuccessUseQuotaWdlConsumedValue() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quotaConsumed to 2000
    final Integer quotaConsumedAbovePipelineMinQuota = 2000;
    final Map<String, String> quotaOutputs =
        new HashMap<>(Map.of("quotaConsumed", quotaConsumedAbovePipelineMinQuota.toString()));
    flightContext.getWorkingMap().put(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // before running make sure quota consumed for user is 0
    UserQuotaTotals userQuotaTotals =
        quotasService.getQuotaTotalsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(0, userQuotaTotals.totalConsumed());

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // after running make sure quota for user is 3000 based on the quota consumed
    // from the job running
    userQuotaTotals =
        quotasService.getQuotaTotalsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(quotaConsumedAbovePipelineMinQuota, userQuotaTotals.totalConsumed());

    // make sure the raw quota assumed was saved for the pipeline
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_1_ID)
            .orElseThrow();
    assertEquals(quotaConsumedAbovePipelineMinQuota, updatedPipelineRun.getRawQuotaConsumed());

    // assert raw quota consumed are correctly stored in the working map
    assertEquals(
        quotaConsumedAbovePipelineMinQuota,
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.RAW_QUOTA_CONSUMED, Integer.class));

    // assert that the quota consumption records are stored in the working map
    List<QuotaConsumptionRecord> records =
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, List.class);

    assertEquals(1, records.size());

    QuotaConsumptionRecord record = records.get(0);
    assertEquals(2000, record.amountConsumed());
  }

  @Test
  void doStepOverQuotaLimit() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quota consumed in working map to above the pipeline quota limit
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "11000"));
    flightContext.getWorkingMap().put(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());

    // make sure we returned the correct error message for the insufficient quota error
    assertEquals(
        "User quota exceeded for pipeline array_imputation. User quota limit: 2500, "
            + "Quota consumed before this run: 0, Quota consumed for this run: 11000. "
            + "Please visit https://services.terra.bio/pipelines/quotas to purchase additional quota.",
        result.getException().get().getMessage());

    // make sure the raw quota assumed was saved for the pipeline
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_1_ID)
            .orElseThrow();
    assertEquals(11000, updatedPipelineRun.getRawQuotaConsumed());
  }

  @Test
  void doStepOverMaximumQuota() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quota consumed in working map to above the pipeline quota limit
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "50000"));
    flightContext.getWorkingMap().put(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());

    // make sure we returned the correct error message
    assertEquals(
        "The number of submitted samples for this job (50000) exceeds the maximum allowed (20000). Please try again with a smaller input size.",
        result.getException().get().getMessage());
  }

  @Test
  void doStepMissingQuotaField() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // missing quota consumed in working map
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of());
    flightContext.getWorkingMap().put(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "An error occurred while running the job. Quota consumption information not found. Please contact support for help.",
        result.getException().get().getMessage());
  }

  @Test
  void doStepNegativeQuota() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quota consumed in working map to a negative value
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "-5"));
    flightContext.getWorkingMap().put(WdlBasedPipelineJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "An error occurred while running the job. No valid input files found; please check your input. Please contact support for help.",
        result.getException().get().getMessage());
  }

  @Test
  void doStepSuccessConsumesAcrossMultipleAllocations() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // add a second allocation so that consuming 3000 spans both of the user's allocations
    QuotaAllocation paidAllocation =
        quotasService.allocateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_1_ID,
            PipelinesEnum.ARRAY_IMPUTATION,
            QuotaAllocationSourceEnum.EXTERNAL_PAID,
            1000,
            null);

    QuotaAllocation defaultAllocation =
        quotasService
            .getAllQuotaAllocationsForUserAndPipeline(
                TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .stream()
            .filter(a -> a.getQuotaSource() == QuotaAllocationSourceEnum.DEFAULT_FREE)
            .findFirst()
            .orElseThrow();
    putQuotaOutputs("3000");

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // the oldest allocation is consumed first, then the next one, and both are in the receipt
    List<QuotaConsumptionRecord> records =
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, List.class);
    assertEquals(2, records.size());
    assertEquals(defaultAllocation.getId(), records.get(0).quotaAllocationId());
    assertEquals(2500, records.get(0).amountConsumed());
    assertEquals(paidAllocation.getId(), records.get(1).quotaAllocationId());
    assertEquals(500, records.get(1).amountConsumed());

    assertEquals(
        3000,
        quotasService
            .getQuotaTotalsForUserAndPipeline(
                TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .totalConsumed());
  }

  @Test
  void doStepReturnsQuotaExceededMessageWhenQuotaRunsOutDuringConsume() {
    // the totals check passes but the consume runs out of quota, e.g. because a
    // concurrent job consumed it in between
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    QuotasService mockQuotasService = mockQuotasServiceWithQuotaAvailable();
    when(mockQuotasService.consumeQuotaFromAllocationsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 1000))
        .thenThrow(
            new BadRequestException(
                "Insufficient quota available to consume for user testUser and pipeline ARRAY_IMPUTATION. Requested: 1000, available: 0"));
    putQuotaOutputs("1000");

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(mockQuotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // the user sees the same friendly message as when the totals check fails, not the internal one
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "User quota exceeded for pipeline array_imputation. User quota limit: 2500, "
            + "Quota consumed before this run: 0, Quota consumed for this run: 1000. "
            + "Please visit https://services.terra.bio/pipelines/quotas to purchase additional quota.",
        result.getException().get().getMessage());
    assertNothingStoredInWorkingMap();
  }

  @Test
  void doStepRetriesOnTransientDatabaseErrorDuringQuotaConsume() {
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    QuotasService mockQuotasService = mockQuotasServiceWithQuotaAvailable();
    when(mockQuotasService.consumeQuotaFromAllocationsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 1000))
        .thenThrow(new CannotAcquireLockException("could not obtain lock on row"));
    putQuotaOutputs("1000");

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(mockQuotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // the consume transaction was rolled back, so Stairway is asked to retry the step
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
    assertNothingStoredInWorkingMap();
  }

  @Test
  void doStepFailsWithGenericMessageOnDataIntegrityErrorDuringConsume() {
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    QuotasService mockQuotasService = mockQuotasServiceWithQuotaAvailable();
    when(mockQuotasService.consumeQuotaFromAllocationsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 1000))
        .thenThrow(
            new InternalServerErrorException(
                "Data integrity issue: allocation 5 for user testUser pipeline ARRAY_IMPUTATION is ACTIVE but has 0 or negative quota remaining (0). Aborting quota consumption."));
    putQuotaOutputs("1000");

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(mockQuotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // the user sees a generic message that does not expose the internal details
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "An error occurred while running the job. Error updating quota for job. Please contact support for help.",
        result.getException().get().getMessage());
    assertNothingStoredInWorkingMap();
  }

  @Test
  void doStepFailsWhenNoQuotaConsumptionRecordsAreReturned() {
    // this should never happen, but the step must not continue without a receipt
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    QuotasService mockQuotasService = mockQuotasServiceWithQuotaAvailable();
    when(mockQuotasService.consumeQuotaFromAllocationsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 1000))
        .thenReturn(List.of());
    putQuotaOutputs("1000");

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(mockQuotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "An error occurred while running the job. Something went wrong while updating quota for user %s and pipeline array_imputation for flight %s. Please contact support for help."
            .formatted(TestUtils.TEST_USER_1_ID, TestUtils.TEST_NEW_UUID),
        result.getException().get().getMessage());
    assertNothingStoredInWorkingMap();
  }

  @Test
  void undoStepWithoutRecordsInWorkingMapDoesNotRevertAnyQuota() {
    // there is no receipt of quota consumption in working map
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    quotasService.consumeQuotaFromAllocationsForUserAndPipeline(
        TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 600);

    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        600,
        quotasService
            .getQuotaTotalsForUserAndPipeline(
                TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .totalConsumed());
  }

  @Test
  void undoStepRetriesOnTransientDatabaseErrorDuringRevert() {
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    flightContext
        .getWorkingMap()
        .put(
            WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS,
            new ArrayList<>(List.of(new QuotaConsumptionRecord(1L, 100))));
    QuotasService mockQuotasService = mock(QuotasService.class);
    doThrow(new CannotAcquireLockException("could not obtain lock on row"))
        .when(mockQuotasService)
        .revertQuotaConsumption(any());

    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(mockQuotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.undoStep(flightContext);

    // the revert transaction was rolled back, so Stairway is asked to retry the undo
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void undoStepFailsWithGenericMessageWhenRevertFails() {
    // consumption records references an allocation that does not exist
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    flightContext
        .getWorkingMap()
        .put(
            WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS,
            new ArrayList<>(List.of(new QuotaConsumptionRecord(-1L, 100))));

    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.undoStep(flightContext);

    // the user is not shown internal details about the allocation that was not found
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertEquals(
        "An error occurred while running the job. Error while updating quota for job. Please contact support for help.",
        result.getException().get().getMessage());
  }

  @Test
  void undoStepSuccess() {
    // setup user
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // add a second allocation and consume 600 quota across the user's allocations
    quotasService.allocateQuotaForUserAndPipeline(
        TestUtils.TEST_USER_1_ID,
        PipelinesEnum.ARRAY_IMPUTATION,
        QuotaAllocationSourceEnum.EXTERNAL_PAID,
        1000,
        null);

    List<QuotaConsumptionRecord> records =
        quotasService.consumeQuotaFromAllocationsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 600);

    flightContext
        .getWorkingMap()
        .put(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, records);

    assertEquals(
        600,
        quotasService
            .getQuotaTotalsForUserAndPipeline(
                TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .totalConsumed());

    // set up the pipeline run with raw quota consumed of 500
    UUID undoStepTestJobId = TestUtils.TEST_NEW_UUID_2;
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(undoStepTestJobId);
    pipelineRun.setRawQuotaConsumed(500);
    pipelineRunsRepository.save(pipelineRun);

    // run the undo step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.undoStep(flightContext);
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // the consumption records are reverted, so no quota is consumed
    UserQuotaTotals totals =
        quotasService.getQuotaTotalsForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(0, totals.totalConsumed());

    // make sure the raw quota assumed was still saved for the pipeline (we didn't undo this part)
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(undoStepTestJobId, TestUtils.TEST_USER_1_ID)
            .orElseThrow();
    assertEquals(500, updatedPipelineRun.getRawQuotaConsumed());
  }
}
