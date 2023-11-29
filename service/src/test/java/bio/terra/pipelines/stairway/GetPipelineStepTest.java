package bio.terra.pipelines.stairway;

import static java.lang.Boolean.FALSE;
import static java.lang.Boolean.TRUE;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

public class GetPipelineStepTest extends BaseContainerTest {

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
  void getPipeline_success() throws InterruptedException {
    String pipelineId = "imputation";
    FlightMap inputParameters = flightContext.getInputParameters();
    inputParameters.put(GetPipelineFlightMapKeys.PIPELINE_ID, pipelineId);

    var getPipelineStep = new GetPipelineStep(pipelinesService, inputParameters);
    Pipeline pipelineInfoResult = pipelinesService.getPipeline(pipelineId);

    var result = getPipelineStep.doStep(flightContext);

    assertEquals(
        pipelineId,
        flightContext.getInputParameters().get(GetPipelineFlightMapKeys.PIPELINE_ID, String.class));

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(TRUE, flightContext.getWorkingMap().get("updateComplete", Boolean.class));
    assertEquals(
        pipelineInfoResult,
        flightContext
            .getWorkingMap()
            .get(StairwayJobMapKeys.RESPONSE.getKeyName(), Pipeline.class));
  }

  //  @Test
  //  void getPipeline_notFound() throws InterruptedException {
  //    String pipelineId = "non-existent pipeline";
  //    FlightMap inputParameters = flightContext.getInputParameters();
  //    inputParameters.put(GetPipelineFlightMapKeys.PIPELINE_ID, pipelineId);
  //
  //    var getPipelineStep = new GetPipelineStep(pipelinesService, inputParameters);
  //
  //    var result = getPipelineStep.doStep(flightContext);
  //
  //    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
  //    assertEquals(FALSE, flightContext.getWorkingMap().get("updateComplete", Boolean.class));
  //  }

  @Test
  void getPipeline_undo() throws InterruptedException {
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put("updateComplete", FALSE);

    var getPipelineStep = new GetPipelineStep(pipelinesService, flightContext.getInputParameters());

    var result = getPipelineStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
