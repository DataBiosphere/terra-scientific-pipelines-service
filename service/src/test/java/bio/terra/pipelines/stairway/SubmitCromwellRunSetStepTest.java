package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.MethodDetails;
import bio.terra.cbas.model.MethodListResponse;
import bio.terra.cbas.model.MethodVersionDetails;
import bio.terra.cbas.model.RunSetStateResponse;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.pipelines.stairway.imputation.SubmitCromwellRunSetStep;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.*;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class SubmitCromwellRunSetStepTest extends BaseEmbeddedDbTest {
  @Mock private CbasService cbasService;
  @Mock private SamService samService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();
    workingMap.put(RunImputationJobFlightMapKeys.CBAS_URI, "cbasUri");

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    UUID runSetId = UUID.randomUUID();
    MethodListResponse getAllMethodsResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name(
                        flightContext
                            .getInputParameters()
                            .get(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, String.class))
                    .addMethodVersionsItem(
                        new MethodVersionDetails().methodVersionId(UUID.randomUUID())));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.checkHealth(any(), any()))
        .thenReturn(new HealthCheckWorkspaceApps.Result(true, "cbas is healthy"));
    when(cbasService.getAllMethods(any(), any())).thenReturn(getAllMethodsResponse);
    when(cbasService.createRunSet(any(), any(), any()))
        .thenReturn(new RunSetStateResponse().runSetId(runSetId));

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        runSetId,
        flightContext.getWorkingMap().get(RunImputationJobFlightMapKeys.RUN_SET_ID, UUID.class));
  }

  @Test
  void doStepUnhealthyCbas() throws InterruptedException {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.checkHealth(any(), any()))
        .thenReturn(new HealthCheckWorkspaceApps.Result(false, "wds is not healthy"));

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // make sure the appropriate step status was returned
    assertEquals(StepStatus.STEP_RESULT_FAILURE_RETRY, result.getStepStatus());
  }

  @Test
  void doStepNoMatchingMethod() throws InterruptedException {
    // setup
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    MethodListResponse getAllMethodsResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name("random name that doesnt match anything")
                    .addMethodVersionsItem(
                        new MethodVersionDetails().methodVersionId(UUID.randomUUID())));
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(cbasService.checkHealth(any(), any()))
        .thenReturn(new HealthCheckWorkspaceApps.Result(true, "cbas is healthy"));
    when(cbasService.getAllMethods(any(), any())).thenReturn(getAllMethodsResponse);

    // do the step
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService);
    StepResult result = submitCromwellRunSetStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }

  @Test
  void undoStepSuccess() throws InterruptedException {
    SubmitCromwellRunSetStep submitCromwellRunSetStep =
        new SubmitCromwellRunSetStep(cbasService, samService);
    StepResult result = submitCromwellRunSetStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
