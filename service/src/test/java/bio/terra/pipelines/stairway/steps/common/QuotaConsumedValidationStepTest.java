package bio.terra.pipelines.stairway.steps.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.UserQuotasRepository;
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
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class QuotaConsumedValidationStepTest extends BaseEmbeddedDbTest {

  @Autowired private QuotasService quotasService;
  @Autowired private UserQuotasRepository userQuotasRepository;

  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
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
        new QuotaConsumedValidationStep(quotasService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // after running make sure quota for user is 500 from the pipeline min quota
    userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    PipelineQuota pipelineQuota = quotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(pipelineQuota.getMinQuotaConsumed(), userQuota.getQuotaConsumed());

    // assert effective quota consumed is correctly stored in the working map
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
        new QuotaConsumedValidationStep(quotasService);
    StepResult result = quotaConsumedValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // after running make sure quota for user is 3000 based on the quota consumed
    // from the job running
    userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(2000, userQuota.getQuotaConsumed());
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
        new QuotaConsumedValidationStep(quotasService);
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
        new QuotaConsumedValidationStep(quotasService);
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

    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService);
    StepResult result = quotaConsumedValidationStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    // the working map has 600 effective quota consumed, so the user quota consumed should be 600 -
    // 500
    // = 100
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(100, userQuota.getQuotaConsumed());
  }
}
