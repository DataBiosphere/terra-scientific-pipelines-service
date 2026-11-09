package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.testutils.TestUtils.CONTROL_WORKSPACE_ID;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.db.entities.QuotaAllocationDeduction;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.db.repositories.QuotaAllocationDeductionsRepository;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.model.QuotaConsumptionRecord;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.stream.Collectors;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class CompletePipelineRunStepTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Autowired private PipelineOutputsRepository pipelineOutputsRepository;
  @Autowired private QuotasService quotasService;
  @Autowired private QuotaAllocationDeductionsRepository quotaAllocationDeductionsRepository;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final Integer effectiveQuotaConsumed = 500;
  private FlightMap inputParameters;
  private FlightMap workingMap;

  private QuotaAllocation firstAllocation;
  private QuotaAllocation secondAllocation;
  private QuotaAllocation thirdAllocation;

  private void assertQuotaDeductionsRecorded(Long pipelineRunId) {
    List<QuotaAllocationDeduction> deductions =
        quotaAllocationDeductionsRepository.findByPipelineRunId(pipelineRunId);

    // one deduction per consumption record, adding up to the quota consumed for the run
    assertEquals(3, deductions.size());
    assertEquals(500, deductions.stream().mapToInt(QuotaAllocationDeduction::getAmount).sum());

    // each allocation was deducted from by exactly the amount in its consumption record
    Map<Long, Integer> amountByAllocationId =
        deductions.stream()
            .collect(
                Collectors.toMap(
                    QuotaAllocationDeduction::getQuotaAllocationId,
                    QuotaAllocationDeduction::getAmount));
    assertEquals(
        Map.of(
            firstAllocation.getId(), 100,
            secondAllocation.getId(), 200,
            thirdAllocation.getId(), 200),
        amountByAllocationId);

    // every deduction is tied to this run
    assertTrue(deductions.stream().allMatch(d -> pipelineRunId.equals(d.getPipelineRunId())));
  }

  @BeforeEach
  void setup() {
    inputParameters = new FlightMap();
    workingMap = new FlightMap();

    // the deductions table has a foreign key to quota_allocations, so the quota consumption
    // records must reference allocations that exist
    PipelinesEnum pipeline = PipelinesEnum.ARRAY_IMPUTATION;
    firstAllocation =
        quotasService.allocateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, pipeline, QuotaAllocationSourceEnum.DEFAULT_FREE, 100, null);
    secondAllocation =
        quotasService.allocateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_1_ID,
            pipeline,
            QuotaAllocationSourceEnum.DEVELOPMENT_FREE,
            200,
            null);
    thirdAllocation =
        quotasService.allocateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, pipeline, QuotaAllocationSourceEnum.EXTERNAL_PAID, 200, null);

    List<QuotaConsumptionRecord> quotaConsumptionRecords =
        new ArrayList<>(
            List.of(
                new QuotaConsumptionRecord(firstAllocation.getId(), 100),
                new QuotaConsumptionRecord(secondAllocation.getId(), 200),
                new QuotaConsumptionRecord(thirdAllocation.getId(), 200)));

    workingMap.put(
        WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS, TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE);
    workingMap.put(WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS, quotaConsumptionRecords);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccessWithFileSizes() {
    // setup
    workingMap.put(
        WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE,
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_SIZE);
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // write the run to the db
    pipelineRunsRepository.save(
        new PipelineRun(
            testJobId,
            TestUtils.TEST_USER_1_ID,
            TestUtils.TEST_PIPELINE_KEY_1,
            TestUtils.TEST_TOOL_VERSION_1,
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
            TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
            null,
            null,
            CommonPipelineRunStatusEnum.SUCCEEDED,
            TestUtils.TEST_PIPELINE_DESCRIPTION_1,
            null,
            null));

    // do the step
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService, quotasService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // make sure the run was updated with isSuccess and quotaConsumed
    PipelineRun writtenJob =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, inputParams.get(JobMapKeys.USER_ID, String.class))
            .orElseThrow();
    assertEquals(CommonPipelineRunStatusEnum.SUCCEEDED, writtenJob.getStatus());
    assertEquals(effectiveQuotaConsumed, writtenJob.getQuotaConsumed());
    assertTrue(writtenJob.getStatus().isSuccess());

    // verify quota deductions were recorded for the run
    assertQuotaDeductionsRecorded(writtenJob.getId());

    List<PipelineOutput> pipelineOutputList =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(writtenJob.getId());
    assertEquals(2, pipelineOutputList.size());

    // assert file output values were written correctly to the database
    PipelineOutput fileOutput =
        pipelineOutputList.stream()
            .filter(output -> output.getOutputName().equals("testFileOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals(
        "gs://fc-secure-%s/testFileOutputValue".formatted(CONTROL_WORKSPACE_ID),
        fileOutput.getOutputValue());
    assertEquals(256, fileOutput.getFileSizeBytes());

    // assert string output values were written correctly to the database
    PipelineOutput stringOutput =
        pipelineOutputList.stream()
            .filter(output -> output.getOutputName().equals("testStringOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals("testStringOutputValue", stringOutput.getOutputValue());
    assertNull(stringOutput.getFileSizeBytes());
  }

  @Test
  void doStepSuccessWithoutFileSizes() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // write the run to the db
    pipelineRunsRepository.save(
        new PipelineRun(
            testJobId,
            TestUtils.TEST_USER_1_ID,
            TestUtils.TEST_PIPELINE_KEY_1,
            TestUtils.TEST_TOOL_VERSION_1,
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
            TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
            null,
            null,
            CommonPipelineRunStatusEnum.SUCCEEDED,
            TestUtils.TEST_PIPELINE_DESCRIPTION_1,
            null,
            null));

    // do the step
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService, quotasService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // make sure the run was updated with isSuccess and quotaConsumed
    PipelineRun writtenJob =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, inputParams.get(JobMapKeys.USER_ID, String.class))
            .orElseThrow();
    assertEquals(CommonPipelineRunStatusEnum.SUCCEEDED, writtenJob.getStatus());
    assertEquals(effectiveQuotaConsumed, writtenJob.getQuotaConsumed());
    assertTrue(writtenJob.getStatus().isSuccess());

    // verify quota deductions were recorded for the run
    assertQuotaDeductionsRecorded(writtenJob.getId());

    List<PipelineOutput> pipelineOutputList =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(writtenJob.getId());
    assertEquals(2, pipelineOutputList.size());

    // assert file output values were written correctly to the database
    PipelineOutput fileOutput =
        pipelineOutputList.stream()
            .filter(output -> output.getOutputName().equals("testFileOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals(
        "gs://fc-secure-%s/testFileOutputValue".formatted(CONTROL_WORKSPACE_ID),
        fileOutput.getOutputValue());
    assertNull(fileOutput.getFileSizeBytes());

    // assert string output values were written correctly to the database
    PipelineOutput stringOutput =
        pipelineOutputList.stream()
            .filter(output -> output.getOutputName().equals("testStringOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals("testStringOutputValue", stringOutput.getOutputValue());
    assertNull(stringOutput.getFileSizeBytes());
  }

  @Test
  void doStepSucceedsWhenRecordingQuotaDeductionsFails() {
    // the records reference an allocation that does not exist, so writing the
    // deductions fails with a foreign key violation
    workingMap.put(
        WdlBasedPipelineJobMapKeys.QUOTA_CONSUMPTION_RECORDS,
        new ArrayList<>(List.of(new QuotaConsumptionRecord(-1L, 500))));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // write the run to the db
    pipelineRunsRepository.save(
        new PipelineRun(
            testJobId,
            TestUtils.TEST_USER_1_ID,
            TestUtils.TEST_PIPELINE_KEY_1,
            TestUtils.TEST_TOOL_VERSION_1,
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
            TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
            null,
            null,
            CommonPipelineRunStatusEnum.SUCCEEDED,
            TestUtils.TEST_PIPELINE_DESCRIPTION_1,
            null,
            null));

    // do the step
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService, quotasService);
    var result = writeJobStep.doStep(flightContext);

    // failure writing the deductions rows must not fail the step for a run that succeeded
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    PipelineRun writtenJob =
        pipelineRunsRepository
            .findByJobIdAndUserId(
                testJobId, flightContext.getInputParameters().get(JobMapKeys.USER_ID, String.class))
            .orElseThrow();
    assertEquals(CommonPipelineRunStatusEnum.SUCCEEDED, writtenJob.getStatus());
    assertEquals(effectiveQuotaConsumed, writtenJob.getQuotaConsumed());

    // the deductions were not recorded since the allocation did not exist
    assertTrue(
        quotaAllocationDeductionsRepository.findByPipelineRunId(writtenJob.getId()).isEmpty());
  }

  @Test
  void undoStepSuccess() {
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService, quotasService);
    var result = writeJobStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
