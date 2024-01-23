package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.Stairway;
import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
import java.util.concurrent.TimeoutException;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class RetryUtilsTest extends BaseContainerTest {

  // defaults give enough time for one retry to complete
  static final Duration TEST_RETRY_TOTAL_DURATION = Duration.ofSeconds(2);
  static final Duration TEST_RETRY_SLEEP_DURATION = Duration.ofSeconds(1);
  static final double TEST_RETRY_FACTOR_INCREASE = 1;
  static final Duration TEST_RETRY_SLEEP_DURATION_MAX = Duration.ofSeconds(1);

  final UUID testFlightId = TestUtils.TEST_NEW_UUID;
  final String testFlightIdString = testFlightId.toString();

  final FlightState runningFlightState =
      StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, testFlightId);
  final FlightState successFlightState =
      StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, testFlightId);

  @Mock private Stairway mockStairway;

  @Test
  void getWithRetryOnException_successOnFirstTry() throws Exception {
    when(mockStairway.getFlightState(testFlightIdString)).thenReturn(successFlightState);

    FlightState flightStateResult =
        RetryUtils.getWithRetryOnException(
            () -> mockStairway.getFlightState(testFlightIdString),
            TEST_RETRY_TOTAL_DURATION,
            TEST_RETRY_SLEEP_DURATION,
            TEST_RETRY_FACTOR_INCREASE,
            TEST_RETRY_SLEEP_DURATION_MAX,
            null);

    assertEquals(FlightStatus.SUCCESS, flightStateResult.getFlightStatus());
  }

  @Test
  void getWithRetryOnException_successAfterRetry() throws Exception {

    InternalStairwayException expectedException = new InternalStairwayException("test exception");

    when(mockStairway.getFlightState(testFlightIdString))
        .thenThrow(expectedException)
        .thenReturn(successFlightState);

    FlightState flightStateResult =
        RetryUtils.getWithRetryOnException(
            () -> mockStairway.getFlightState(testFlightIdString),
            TEST_RETRY_TOTAL_DURATION,
            TEST_RETRY_SLEEP_DURATION,
            TEST_RETRY_FACTOR_INCREASE,
            TEST_RETRY_SLEEP_DURATION_MAX,
            null);

    assertEquals(FlightStatus.SUCCESS, flightStateResult.getFlightStatus());
  }

  @Test
  void getWithRetryOnException_retryableIsRetried() throws Exception {
    // if the exception is not in the list of retryable exceptions, it should be thrown
    List<Class<? extends Exception>> retryableExceptions = new ArrayList<>();
    retryableExceptions.add(InvalidResultStateException.class);

    InvalidResultStateException retryableException =
        new InvalidResultStateException("test exception");

    when(mockStairway.getFlightState(testFlightIdString))
        .thenThrow(retryableException)
        .thenReturn(successFlightState);

    FlightState flightStateResult =
        RetryUtils.getWithRetryOnException(
            () -> mockStairway.getFlightState(testFlightIdString),
            TEST_RETRY_TOTAL_DURATION,
            TEST_RETRY_SLEEP_DURATION,
            TEST_RETRY_FACTOR_INCREASE,
            TEST_RETRY_SLEEP_DURATION_MAX,
            retryableExceptions);

    assertEquals(FlightStatus.SUCCESS, flightStateResult.getFlightStatus());
  }

  @Test
  void getWithRetryOnException_nonRetryableThrows() throws Exception {
    // if the exception is not in the list of retryable exceptions, it should be thrown
    List<Class<? extends Exception>> retryableExceptions = new ArrayList<>();
    retryableExceptions.add(InvalidResultStateException.class);

    InternalStairwayException nonRetryableException =
        new InternalStairwayException("test exception");

    when(mockStairway.getFlightState(testFlightIdString))
        .thenThrow(nonRetryableException)
        .thenReturn(successFlightState);

    assertThrows(
        InternalStairwayException.class,
        () ->
            RetryUtils.getWithRetryOnException(
                () -> mockStairway.getFlightState(testFlightIdString),
                TEST_RETRY_TOTAL_DURATION,
                TEST_RETRY_SLEEP_DURATION,
                TEST_RETRY_FACTOR_INCREASE,
                TEST_RETRY_SLEEP_DURATION_MAX,
                retryableExceptions));
  }

  @Test
  void getWithRetryOnException_failAfterRetryTimeout() throws Exception {
    // make the retries take longer than the total allowed duration
    Duration shortTotalDuration = Duration.ofSeconds(1);
    Duration longRetrySleepDuration = Duration.ofSeconds(2);

    InternalStairwayException expectedException = new InternalStairwayException("test exception");

    // this should fail once, wait a second and retry, fail again, and run out of time
    when(mockStairway.getFlightState(testFlightIdString))
        .thenThrow(expectedException)
        .thenThrow(expectedException)
        .thenReturn(successFlightState);

    assertThrows(
        InternalStairwayException.class,
        () ->
            RetryUtils.getWithRetryOnException(
                () -> mockStairway.getFlightState(testFlightIdString),
                shortTotalDuration,
                longRetrySleepDuration,
                TEST_RETRY_FACTOR_INCREASE,
                TEST_RETRY_SLEEP_DURATION_MAX,
                null));
  }

  // TODO test with list of allowable exceptions

  @Test
  void getWithRetry_successAfterOneRetry() throws Exception {

    when(mockStairway.getFlightState(testFlightIdString))
        .thenReturn(runningFlightState)
        .thenReturn(successFlightState);

    FlightState flightStateResult =
        RetryUtils.getWithRetry(
            FlightUtils::flightComplete,
            () -> mockStairway.getFlightState(testFlightIdString),
            TEST_RETRY_TOTAL_DURATION,
            TEST_RETRY_SLEEP_DURATION,
            TEST_RETRY_FACTOR_INCREASE,
            TEST_RETRY_SLEEP_DURATION_MAX);

    assertEquals(FlightStatus.SUCCESS, flightStateResult.getFlightStatus());
  }

  @Test
  void getWithRetry_timeoutAfterRetry() throws Exception {
    // make the retries take longer than the total allowed duration
    Duration shortTotalDuration = Duration.ofSeconds(1);
    Duration longRetrySleepDuration = Duration.ofSeconds(2);

    // this should check once, retry after 2 seconds, check again, and run out of time
    when(mockStairway.getFlightState(testFlightIdString))
        .thenReturn(runningFlightState)
        .thenReturn(runningFlightState)
        .thenReturn(successFlightState);

    assertThrows(
        TimeoutException.class,
        () ->
            RetryUtils.getWithRetry(
                FlightUtils::flightComplete,
                () -> mockStairway.getFlightState(testFlightIdString),
                shortTotalDuration,
                longRetrySleepDuration,
                TEST_RETRY_FACTOR_INCREASE,
                TEST_RETRY_SLEEP_DURATION_MAX));
  }
}
