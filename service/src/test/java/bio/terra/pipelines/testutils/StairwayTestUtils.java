package bio.terra.pipelines.testutils;

import static bio.terra.stairway.FlightStatus.*;
import static org.testcontainers.shaded.org.awaitility.Awaitility.await;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.DatabaseOperationException;
import bio.terra.stairway.exception.DuplicateFlightIdException;
import bio.terra.stairway.exception.StairwayExecutionException;
import java.time.Instant;
import java.util.HashMap;
import java.util.List;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** Test utilities for working with Stairway. */
public class StairwayTestUtils {
  private static final Logger logger = LoggerFactory.getLogger(StairwayTestUtils.class);

  private StairwayTestUtils() {}

  public static final Instant TIME_SUBMITTED_1 = Instant.parse("2024-01-01T00:00:00.00Z");
  public static final Instant TIME_SUBMITTED_2 = Instant.parse("2024-01-02T01:00:00.00Z");
  public static final Instant TIME_COMPLETED_1 = Instant.parse("2024-01-01T00:30:00.00Z");
  public static final Instant TIME_COMPLETED_2 = Instant.parse("2024-01-02T01:30:00.00Z");
  public static final FlightMap CREATE_JOB_INPUT_PARAMS =
      StairwayTestUtils.constructCreateJobInputs(
          TestUtils.TEST_PIPELINE_1_ENUM,
          TestUtils.TEST_PIPELINE_ID_1,
          TestUtils.TEST_USER_ID_1,
          TestUtils.TEST_PIPELINE_INPUTS,
          TestUtils.CONTROL_WORKSPACE_ID,
          TestUtils.TEST_WDL_METHOD_NAME_1,
              TestUtils.TEST_RESULT_URL);
  public static final FlightMap EMPTY_WORKING_MAP = new FlightMap();
  public static final String TEST_DESCRIPTION = "Test Job Description";

  public static final FlightState FLIGHT_STATE_DONE_SUCCESS_1 =
      StairwayTestUtils.constructFlightStateWithStatusAndId(
          FlightStatus.SUCCESS,
          TestUtils.TEST_NEW_UUID,
          CREATE_JOB_INPUT_PARAMS,
          EMPTY_WORKING_MAP,
          TIME_SUBMITTED_1,
          TIME_COMPLETED_1);
  public static final FlightState FLIGHT_STATE_DONE_SUCCESS_2 =
      StairwayTestUtils.constructFlightStateWithStatusAndId(
          FlightStatus.SUCCESS,
          TestUtils.TEST_NEW_UUID_2,
          CREATE_JOB_INPUT_PARAMS,
          EMPTY_WORKING_MAP,
          TIME_SUBMITTED_2,
          TIME_COMPLETED_2);

  public static final String PAGE_TOKEN = "foo";

  public static final EnumeratedJob ENUMERATED_JOB_DONE_SUCCESS_1 =
      new EnumeratedJob().flightState(FLIGHT_STATE_DONE_SUCCESS_1);
  public static final EnumeratedJob ENUMERATED_JOB_DONE_SUCCESS_2 =
      new EnumeratedJob().flightState(FLIGHT_STATE_DONE_SUCCESS_2);

  public static final EnumeratedJobs ENUMERATED_JOBS =
      new EnumeratedJobs()
          .results(List.of(ENUMERATED_JOB_DONE_SUCCESS_1, ENUMERATED_JOB_DONE_SUCCESS_2))
          .totalResults(2)
          .pageToken(PAGE_TOKEN);

  /**
   * Submits the flight and block until Stairway completes it by polling regularly until the timeout
   * is reached.
   */
  public static FlightState blockUntilFlightCompletes(
      Stairway stairway,
      Class<? extends Flight> flightClass,
      UUID flightId,
      FlightMap inputParameters,
      Long timeoutInSeconds,
      FlightDebugInfo debugInfo)
      throws DatabaseOperationException, StairwayExecutionException, InterruptedException,
          DuplicateFlightIdException {
    stairway.submitWithDebugInfo(
        flightId.toString(), flightClass, inputParameters, /* shouldQueue= */ false, debugInfo);
    return pollUntilComplete(flightId, stairway, timeoutInSeconds);
  }

  /**
   * Polls stairway until the flight for {@code flightId} completes, or this has polled until
   * {@timeoutInSeconds} seconds have elapsed.
   */
  public static FlightState pollUntilComplete(
      UUID flightId, Stairway stairway, Long timeoutInSeconds)
      throws InterruptedException, DatabaseOperationException {
    String flightIdString = flightId.toString(); // Stairway expects a String flightId
    await()
        .atMost(timeoutInSeconds, TimeUnit.SECONDS)
        .until(() -> !stairway.getFlightState(flightIdString).isActive());
    FlightState flightState = stairway.getFlightState(flightIdString);
    if (!flightState.isActive()) {
      return stairway.getFlightState(flightIdString);
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
      PipelinesEnum pipelineName,
      Long pipelineId,
      String userId,
      Object pipelineInputs,
      UUID controlWorkspaceId,
      String wdlMethodName,
      String resultPath) {
    FlightMap inputParameters = new FlightMap();
    return constructCreateJobInputs(
        inputParameters,
        pipelineName,
        pipelineId,
        userId,
        pipelineInputs,
        controlWorkspaceId,
        wdlMethodName,
            resultPath);
  }

  public static FlightMap constructCreateJobInputs(
      FlightMap inputParameters,
      PipelinesEnum pipelineName,
      Long pipelineId,
      String userId,
      Object pipelineInputs,
      UUID controlWorkspaceId,
      String wdlMethodName,
      String resultPath) {
    inputParameters.put(JobMapKeys.USER_ID.getKeyName(), userId);
    inputParameters.put(JobMapKeys.PIPELINE_NAME.getKeyName(), pipelineName);
    inputParameters.put(JobMapKeys.DESCRIPTION.getKeyName(), TEST_DESCRIPTION);
    inputParameters.put(JobMapKeys.RESULT_PATH.getKeyName(), resultPath);
    inputParameters.put(RunImputationJobFlightMapKeys.PIPELINE_ID, pipelineId);
    inputParameters.put(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs);
    inputParameters.put(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID, controlWorkspaceId);
    inputParameters.put(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, wdlMethodName);

    return inputParameters;
  }

  public static FlightMap constructCreateJobInputs(FlightMap inputParameters) {
    return constructCreateJobInputs(
        inputParameters,
        PipelinesEnum.IMPUTATION_MINIMAC4,
        TestUtils.TEST_PIPELINE_ID_1,
        TestUtils.TEST_USER_ID_1,
        new HashMap<>(),
        TestUtils.CONTROL_WORKSPACE_ID,
        TestUtils.TEST_WDL_METHOD_NAME_1,
            TestUtils.TEST_RESULT_URL);
  }

  /* Construct a FlightState with the given status and id. resultMap and inputParameters will be empty, and timeSubmitted and timeCompleted will be ~now. */
  public static FlightState constructFlightStateWithStatusAndId(
      FlightStatus flightStatus, UUID flightId) {
    FlightMap resultMap = new FlightMap();
    FlightMap inputParameters = constructCreateJobInputs(new FlightMap());
    Instant timeSubmitted = Instant.now();
    Instant timeCompleted = Instant.now();
    return constructFlightStateWithStatusAndId(
        flightStatus, flightId, inputParameters, resultMap, timeSubmitted, timeCompleted);
  }

  /* Construct a FlightState with the given status, id, resultMap, and inputParameters. timeSubmitted and timeCompleted will be ~now. */
  public static FlightState constructFlightStateWithStatusAndId(
      FlightStatus flightStatus, UUID flightId, FlightMap inputParameters, FlightMap resultMap) {
    Instant timeSubmitted = Instant.now();
    Instant timeCompleted = Instant.now();
    return constructFlightStateWithStatusAndId(
        flightStatus, flightId, inputParameters, resultMap, timeSubmitted, timeCompleted);
  }

  public static FlightState constructFlightStateWithStatusAndId(
      FlightStatus flightStatus,
      UUID flightId,
      FlightMap inputParameters,
      FlightMap resultMap,
      Instant submittedTime,
      Instant completedTime) {
    FlightState flightState = new FlightState();
    flightState.setFlightId(flightId.toString());

    flightState.setInputParameters(inputParameters);
    flightState.setResultMap(resultMap);

    flightState.setFlightStatus(flightStatus);
    flightState.setSubmitted(submittedTime);
    if (flightStatus == SUCCESS || flightStatus == ERROR || flightStatus == FATAL) {
      flightState.setCompleted(completedTime);
    }
    return flightState;
  }
}
