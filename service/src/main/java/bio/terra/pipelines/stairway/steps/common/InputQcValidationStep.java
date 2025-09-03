package bio.terra.pipelines.stairway.steps.common;

import bio.terra.common.exception.BadRequestException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.Map;
import java.util.Objects;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step checks that the input qc for the flight is at least the min_quota_consumed for the
 * pipeline being run. Once that is evaluated it then checks that the quota consumed for this run
 * does not cause the user to exceed their quota limit. If everything passes then this step writes
 * the effective quota consumed for this run to the working map.
 *
 * <p>This step expects raw_quota_consumed to be provided in the working map
 */
public class InputQcValidationStep implements Step {

  private final Logger logger = LoggerFactory.getLogger(InputQcValidationStep.class);

  public InputQcValidationStep() {}

  @Override
  @SuppressWarnings(
      "java:S2259") // suppress warning for possible NPE when calling inputQcOutputs.get(),
  //  since we do validate that inputQcOutputs is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) {

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.INPUT_QC_OUTPUTS);

    Map<String, ?> inputQcOutputs =
        workingMap.get(ImputationJobMapKeys.INPUT_QC_OUTPUTS, Map.class);

    // extract passes_qc
    boolean passesQc = (boolean) inputQcOutputs.get("passesQc");
    if (!passesQc) {
      // extract error messages
      String qcMessages = (String) Objects.requireNonNull(inputQcOutputs.get("qcMessages"));
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new BadRequestException("Input failed QC: " + qcMessages));
    }
    logger.info("Input passed QC");
    return StepResult.getStepResultSuccess();
  }

  // Nothing to undo
  @SuppressWarnings("java:S2259") // suppress warning for possible NPE when calling workingMap.get,
  // these values should exist if this step succeeded and is getting undone due to a
  // future step failing.
  @Override
  public StepResult undoStep(FlightContext flightContext) {
    return StepResult.getStepResultSuccess();
  }
}
