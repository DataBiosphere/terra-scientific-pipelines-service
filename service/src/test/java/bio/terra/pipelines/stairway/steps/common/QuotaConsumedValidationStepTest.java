package bio.terra.pipelines.stairway.steps.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.PipelineQuotasRepository;
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
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class QuotaConsumedValidationStepTest extends BaseEmbeddedDbTest {

  @Autowired private QuotasService quotasService;
  @Autowired private UserQuotasRepository userQuotasRepository;
  @Autowired private PipelineQuotasRepository pipelineQuotasRepository;

  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();
    workingMap.put(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, 30);
    workingMap.put(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, 40);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccessUsePipelineQuotaFloor() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

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
  }

  @Test
  void doStepSuccessUseQuotaWdlConsumedValue() {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    flightContext.getWorkingMap().put(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, 2000);

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
    // assert effective quota consumed is correctly stored in the working map
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
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(ImputationJobMapKeys.RAW_QUOTA_CONSUMED, 11000);

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
    // create a user quotas row with 50 quota consumed
    userQuotasRepository.save(
        new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TestUtils.TEST_USER_ID_1, 10000, 50));

    // make sure quota consumed is updated properly for undoStep
    QuotaConsumedValidationStep quotaConsumedValidationStep =
        new QuotaConsumedValidationStep(quotasService);
    StepResult result = quotaConsumedValidationStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    // the working map has 40 effective quota consumed, so the user quota consumed should be 50 - 40
    // = 10
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(10, userQuota.getQuotaConsumed());
  }
}
