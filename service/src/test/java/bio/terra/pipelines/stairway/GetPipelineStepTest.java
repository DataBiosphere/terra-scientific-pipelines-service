package bio.terra.pipelines.stairway;

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
    //            var inputParameters = new FlightMap();
    //            inputParameters.put(
    //                    WorkspaceFlightMapKeys.ControlledResourceKeys.DESTINATION_WORKSPACE_ID,
    //                    destinationWorkspaceId);
    //            inputParameters.put(
    //                    WorkspaceFlightMapKeys.ControlledResourceKeys.PREFIXES_TO_CLONE,
    // clonePrefixes);

    var workingMap = new FlightMap();

    //            when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void getPipeline_success() throws InterruptedException {
    var getPipelineStep = new GetPipelineStep(pipelinesService);
    Pipeline pipelineInfoResult = pipelinesService.getPipeline("imputation");

    var result = getPipelineStep.doStep(flightContext);

    assertEquals(result.getStepStatus(), StepStatus.STEP_RESULT_SUCCESS);
    assertEquals(TRUE, flightContext.getWorkingMap().get("updateComplete", Boolean.class));
    assertEquals(
        pipelineInfoResult,
        flightContext
            .getWorkingMap()
            .get(StairwayJobMapKeys.RESPONSE.getKeyName(), Pipeline.class));
  }
}
