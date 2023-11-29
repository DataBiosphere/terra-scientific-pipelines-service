package bio.terra.pipelines.testutils;

import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightDebugInfo;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.Stairway;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import bio.terra.stairway.exception.DatabaseOperationException;
import bio.terra.stairway.exception.DuplicateFlightIdException;
import bio.terra.stairway.exception.StairwayExecutionException;
import java.time.Duration;
import java.time.Instant;
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
      Duration timeout,
      FlightDebugInfo debugInfo)
      throws DatabaseOperationException, StairwayExecutionException, InterruptedException,
          DuplicateFlightIdException {
    String flightId = stairway.createFlightId();

    stairway.submitWithDebugInfo(
        flightId, flightClass, inputParameters, /* shouldQueue= */ false, debugInfo);
    return pollUntilComplete(flightId, stairway, Duration.ofSeconds(1), timeout);
  }

  /**
   * Polls stairway until the flight for {@code flightId} completes, or this has polled {@code
   * numPolls} times every {@code pollInterval}.
   */
  public static FlightState pollUntilComplete(
      String flightId, Stairway stairway, Duration pollInterval, Duration timeout)
      throws InterruptedException, DatabaseOperationException {
    for (Instant deadline = Instant.now().plus(timeout);
        Instant.now().isBefore(deadline);
        Thread.sleep(pollInterval.toMillis())) {
      FlightState flightState = stairway.getFlightState(flightId);
      if (!flightState.isActive()) {
        return flightState;
      }
    }
    throw new InterruptedException(
        String.format("Flight [%s] did not complete in the allowed wait time.", flightId));
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
}
