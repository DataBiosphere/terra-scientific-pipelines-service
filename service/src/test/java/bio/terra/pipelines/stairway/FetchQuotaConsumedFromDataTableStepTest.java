package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class FetchQuotaConsumedFromDataTableStepTest extends BaseEmbeddedDbTest {

  @Mock RawlsService rawlsService;
  @Mock SamService samService;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("thisToken");
  }

  @Test
  void doStepSuccess() throws RawlsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    // outputs to match the test output definitions
    int quotaConsumed = 15;
    Map<String, Object> entityAttributes = new HashMap<>(Map.of("quota_consumed", quotaConsumed));
    Entity entity = new Entity().attributes(entityAttributes);

    when(rawlsService.getDataTableEntity(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenReturn(entity);

    FetchQuotaConsumedFromDataTableStep fetchQuotaConsumedFromDataTableStep =
        new FetchQuotaConsumedFromDataTableStep(rawlsService, samService);
    StepResult result = fetchQuotaConsumedFromDataTableStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    assertEquals(
        quotaConsumed,
        flightContext.getWorkingMap().get(ImputationJobMapKeys.QUOTA_CONSUMED, Integer.class));
  }

  @Test
  void doStepRawlsFailureRetry() throws RawlsServiceException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    when(rawlsService.getDataTableEntity(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenThrow(new RawlsServiceApiException("Rawls Service Api Exception"));

    FetchQuotaConsumedFromDataTableStep fetchQuotaConsumedFromDataTableStep =
        new FetchQuotaConsumedFromDataTableStep(rawlsService, samService);
    StepResult result = fetchQuotaConsumedFromDataTableStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void doStepOutputsFailureNoRetry() throws InternalServerErrorException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    when(flightContext.getFlightId()).thenReturn(TestUtils.TEST_NEW_UUID.toString());

    // try with no quota_consumed attribute
    Map<String, Object> entityAttributes = new HashMap<>();
    Entity entity = new Entity().attributes(entityAttributes);

    when(rawlsService.getDataTableEntity(
            "thisToken",
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            PipelinesEnum.ARRAY_IMPUTATION.getValue(),
            TestUtils.TEST_NEW_UUID.toString()))
        .thenReturn(entity);

    FetchQuotaConsumedFromDataTableStep fetchQuotaConsumedFromDataTableStep =
        new FetchQuotaConsumedFromDataTableStep(rawlsService, samService);
    StepResult result = fetchQuotaConsumedFromDataTableStep.doStep(flightContext);
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());

    // try with quota_consumed attribute as 0
    entityAttributes.put("quota_consumed", 0);
    result = fetchQuotaConsumedFromDataTableStep.doStep(flightContext);
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() {
    FetchQuotaConsumedFromDataTableStep fetchQuotaConsumedFromDataTableStep =
        new FetchQuotaConsumedFromDataTableStep(rawlsService, samService);
    StepResult result = fetchQuotaConsumedFromDataTableStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
