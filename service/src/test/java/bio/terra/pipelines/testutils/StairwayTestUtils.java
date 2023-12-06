package bio.terra.pipelines.testutils;

import static org.testcontainers.shaded.org.awaitility.Awaitility.await;

import bio.terra.pipelines.stairway.GetPipelineFlightMapKeys;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.DatabaseOperationException;
import bio.terra.stairway.exception.DuplicateFlightIdException;
import bio.terra.stairway.exception.StairwayExecutionException;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** Test utilities for working with Stairway. */
public class StairwayTestUtils {
  private static final Logger logger = LoggerFactory.getLogger(StairwayTestUtils.class);

  private StairwayTestUtils() {}

  /**
   * Submits the flight and block until Stairway completes it by polling regularly until the timeout
   * is reached.
   */
  public static FlightState blockUntilFlightCompletes(
      Stairway stairway,
      Class<? extends Flight> flightClass,
      FlightMap inputParameters,
      Long timeoutInSeconds,
      FlightDebugInfo debugInfo)
      throws DatabaseOperationException, StairwayExecutionException, InterruptedException,
          DuplicateFlightIdException {
    String flightId = stairway.createFlightId();

    stairway.submitWithDebugInfo(
        flightId, flightClass, inputParameters, /* shouldQueue= */ false, debugInfo);
    return pollUntilComplete(flightId, stairway, timeoutInSeconds);
  }

  /**
   * Polls stairway until the flight for {@code flightId} completes, or this has polled until
   * {@timeoutInSeconds} seconds have elapsed.
   */
  public static FlightState pollUntilComplete(
      String flightId, Stairway stairway, Long timeoutInSeconds)
      throws InterruptedException, DatabaseOperationException {
    await()
        .atMost(timeoutInSeconds, TimeUnit.SECONDS)
        .until(() -> !stairway.getFlightState(flightId).isActive());
    FlightState flightState = stairway.getFlightState(flightId);
    if (!flightState.isActive()) {
      return stairway.getFlightState(flightId);
    } else {
      throw new InterruptedException(
          String.format("Flight [%s] did not complete in the allowed wait time.", flightId));
    }
  }

  /**
   * A {@link Step} that always fatally errors on {@link Step#doStep(FlightContext)}. Undo is ok.
   */
  public static class ErrorDoStep implements Step {
    @Override
    public StepResult doStep(FlightContext flightContext) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL);
    }

    @Override
    public StepResult undoStep(FlightContext flightContext) {
      return StepResult.getStepResultSuccess();
    }
  }

  public static FlightMap constructGetPipelineInputs(String pipelineId) {
    FlightMap inputParameters = new FlightMap();
    return constructGetPipelineInputs(inputParameters, pipelineId);
  }

  public static FlightMap constructGetPipelineInputs(FlightMap inputParameters, String pipelineId) {
    inputParameters.put(GetPipelineFlightMapKeys.PIPELINE_ID, pipelineId);
    return inputParameters;
  }

  public static FlightState constructFlightStateWithStatus(
      FlightStatus flightStatus, FlightMap resultMap) {
    FlightState flightState = new FlightState();
    flightState.setFlightId("testFlightId");

    flightState.setResultMap(resultMap);

    flightState.setFlightStatus(flightStatus);
    return flightState;
  }

  public static FlightState constructFlightStateWithStatus(FlightStatus flightStatus) {
    FlightMap resultMap = new FlightMap();
    return constructFlightStateWithStatus(flightStatus, resultMap);
  }
}
