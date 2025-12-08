package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.testutils.TestUtils.createNewPipelineRunWithJobId;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.db.repositories.UserQuotasRepository;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class QuotaConsumedValidationStepTest extends BaseEmbeddedDbTest {

  @Autowired private QuotasService quotasService;
  @Autowired private UserQuotasRepository userQuotasRepository;
  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;

  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

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
  }

  @Test
  void doStepSuccessUsePipelineQuotaFloor() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quotaConsumed to 30 which is below the pipeline min quota of 500
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "30"));
    flightContext.getWorkingMap().put(ImputationJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // before running make sure quota consumed for user is 0
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(0, userQuota.getQuotaConsumed());

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // after running make sure quota for user is 500 from the pipeline min quota
    userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    PipelineQuota pipelineQuota = quotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(pipelineQuota.getMinQuotaConsumed(), userQuota.getQuotaConsumed());

    // make sure the raw quota assumed was saved for the pipeline
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1)
            .orElseThrow();
    assertEquals(30, updatedPipelineRun.getRawQuotaConsumed());

    // assert raw and effective quota consumed are correctly stored in the working map
    assertEquals(
        pipelineQuota.getMinQuotaConsumed(),
        flightContext
            .getWorkingMap()
            .get(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, Integer.class));
    assertEquals(
        30,
        flightContext.getWorkingMap().get(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, Integer.class));
  }

  @Test
  void doStepSuccessUseQuotaWdlConsumedValue() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quotaConsumed to 2000
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "2000"));
    flightContext.getWorkingMap().put(ImputationJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // before running make sure quota consumed for user is 0
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(0, userQuota.getQuotaConsumed());

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // after running make sure quota for user is 3000 based on the quota consumed
    // from the job running
    userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(2000, userQuota.getQuotaConsumed());

    // make sure the raw quota assumed was saved for the pipeline
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1)
            .orElseThrow();
    assertEquals(2000, updatedPipelineRun.getRawQuotaConsumed());

    // assert raw and effective quota consumed are correctly stored in the working map
    assertEquals(
        2000,
        flightContext.getWorkingMap().get(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, Integer.class));
    assertEquals(
        2000,
        flightContext
            .getWorkingMap()
            .get(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, Integer.class));
  }

  @Test
  void doStepOverQuotaLimit() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quota consumed in working map to above the pipeline quota limit
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "11000"));
    flightContext.getWorkingMap().put(ImputationJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());

    // make sure the raw quota assumed was saved for the pipeline
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1)
            .orElseThrow();
    assertEquals(11000, updatedPipelineRun.getRawQuotaConsumed());
  }

  @Test
  void doStepMissingQuotaField() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // missing quota consumed in working map
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of());
    flightContext.getWorkingMap().put(ImputationJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }

  @Test
  void doStepNegativeQuota() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // set quota consumed in working map to a negative value
    final Map<String, String> quotaOutputs = new HashMap<>(Map.of("quotaConsumed", "-5"));
    flightContext.getWorkingMap().put(ImputationJobMapKeys.QUOTA_OUTPUTS, quotaOutputs);

    // do the step
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService, pipelineRunsService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a failure
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    // setup user
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    // create a user quotas row with 1000 quota consumed
    userQuotasRepository.save(
        new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TestUtils.TEST_USER_ID_1, 10000, 600));

    // make sure quota consumed is updated properly for undoStep
    flightContext.getWorkingMap().put(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, 500);

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
    // the user quota consumed was originally 600, and the working map has 500 effective quota
    // consumed,
    // so the user quota consumed should be updated to 600 - 500 = 100
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(100, userQuota.getQuotaConsumed());

    // make sure the raw quota assumed was still saved for the pipeline (we didn't undo this part)
    PipelineRun updatedPipelineRun =
        pipelineRunsRepository
            .findByJobIdAndUserId(undoStepTestJobId, TestUtils.TEST_USER_ID_1)
            .orElseThrow();
    assertEquals(500, updatedPipelineRun.getRawQuotaConsumed());
  }
}
