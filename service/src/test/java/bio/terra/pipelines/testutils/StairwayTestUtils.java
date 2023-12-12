package bio.terra.pipelines.testutils;

import static org.testcontainers.shaded.org.awaitility.Awaitility.await;

import bio.terra.pipelines.stairway.CreateJobFlightMapKeys;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.DatabaseOperationException;
import bio.terra.stairway.exception.DuplicateFlightIdException;
import bio.terra.stairway.exception.StairwayExecutionException;
import java.util.HashMap;
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
      String flightId,
      FlightMap inputParameters,
      Long timeoutInSeconds,
      FlightDebugInfo debugInfo)
      throws DatabaseOperationException, StairwayExecutionException, InterruptedException,
          DuplicateFlightIdException {
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

  public static FlightMap constructCreateJobInputs(
      String pipelineId, String pipelineVersion, String submittingUserId, Object pipelineInputs) {
    FlightMap inputParameters = new FlightMap();
    return constructCreateJobInputs(
        inputParameters, pipelineId, pipelineVersion, submittingUserId, pipelineInputs);
  }

  public static FlightMap constructCreateJobInputs(
      FlightMap inputParameters,
      String pipelineId,
      String pipelineVersion,
      String submittingUserId,
      Object pipelineInputs) {
    inputParameters.put(CreateJobFlightMapKeys.PIPELINE_ID, pipelineId);
    inputParameters.put(CreateJobFlightMapKeys.PIPELINE_VERSION, pipelineVersion);
    inputParameters.put(CreateJobFlightMapKeys.SUBMITTING_USER_ID, submittingUserId);
    inputParameters.put(CreateJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs);

    return inputParameters;
  }

  public static FlightMap constructCreateJobInputs(FlightMap inputParameters) {
    return constructCreateJobInputs(
        inputParameters,
        "testPipelineId",
        "testPipelineVersion",
        "testSubmittingUserId",
        new HashMap<>());
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
