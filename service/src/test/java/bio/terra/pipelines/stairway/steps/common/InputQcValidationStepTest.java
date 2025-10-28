package bio.terra.pipelines.stairway.steps.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class InputQcValidationStepTest extends BaseEmbeddedDbTest {
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void testDoStepPassesQc() {
    Map<String, ?> inputQcOutputs = new HashMap<>(Map.of("passesQc", true, "qcMessages", ""));
    flightContext.getWorkingMap().put(ImputationJobMapKeys.INPUT_QC_OUTPUTS, inputQcOutputs);

    // do the step
    InputQcValidationStep inputQcValidationStep = new InputQcValidationStep();
    StepResult result = inputQcValidationStep.doStep(flightContext);

    // make sure the step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void testDoStepFailsQc() {
    Map<String, ?> inputQcOutputs =
        new HashMap<>(Map.of("passesQc", false, "qcMessages", "File format error."));
    flightContext.getWorkingMap().put(ImputationJobMapKeys.INPUT_QC_OUTPUTS, inputQcOutputs);

    // do the step
    InputQcValidationStep inputQcValidationStep = new InputQcValidationStep();
    StepResult result = inputQcValidationStep.doStep(flightContext);

    // make sure the step failed fatally
    assertEquals(StepStatus.STEP_RESULT_FAILURE_FATAL, result.getStepStatus());
    assertTrue(
        result
            .getException()
            .get()
            .getMessage()
            .contains(
                "User input failed QC: File format error. To troubleshoot, please see documentation at https://broadscientificservices.zendesk.com."));
  }

  @Test
  void testUndoStep() {
    // do the step
    InputQcValidationStep inputQcValidationStep = new InputQcValidationStep();
    StepResult result = inputQcValidationStep.undoStep(flightContext);

    // make sure the undo step was a success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
