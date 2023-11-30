package bio.terra.pipelines.stairway;

import static java.lang.Boolean.FALSE;
import static java.lang.Boolean.TRUE;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class GetPipelineStepTest extends BaseContainerTest {

  @Autowired private PipelinesService pipelinesService;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void getPipeline_doStep_success() throws InterruptedException {
    String pipelineId = "imputation";
    StairwayTestUtils.constructGetPipelineInputs(flightContext.getInputParameters(), pipelineId);

    var getPipelineStep = new GetPipelineStep(pipelinesService);
    Pipeline pipelineInfoResult = pipelinesService.getPipeline(pipelineId);

    var result = getPipelineStep.doStep(flightContext);

    assertEquals(
        pipelineId,
        flightContext.getInputParameters().get(GetPipelineFlightMapKeys.PIPELINE_ID, String.class));

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        TRUE,
        flightContext.getWorkingMap().get(GetPipelineFlightMapKeys.LOOKUP_COMPLETE, Boolean.class));
    assertEquals(
        pipelineInfoResult,
        flightContext
            .getWorkingMap()
            .get(StairwayJobMapKeys.RESPONSE.getKeyName(), Pipeline.class));
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void getPipeline_undoStep_success() throws InterruptedException {
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(GetPipelineFlightMapKeys.LOOKUP_COMPLETE, FALSE);

    var getPipelineStep = new GetPipelineStep(pipelinesService);

    var result = getPipelineStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void getPipeline_undoStep_dismalFailure() throws InterruptedException {
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(GetPipelineFlightMapKeys.LOOKUP_COMPLETE, TRUE);

    var getPipelineStep = new GetPipelineStep(pipelinesService);

    var result = getPipelineStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  }
}
