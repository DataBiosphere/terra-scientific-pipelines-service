package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class JobServiceMockTest extends BaseContainerTest {

  @InjectMocks JobService jobService;
  @Mock private Stairway mockStairway;
  @Mock private StairwayComponent mockStairwayComponent;

  @BeforeEach
  void setup() {
    when(mockStairwayComponent.get()).thenReturn(mockStairway);
  }

  /**
   * Reset the {@link JobService} {@link FlightDebugInfo} after each test so that future submissions
   * aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    jobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void retrieveJobResult_successWithResultClass() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(JobMapKeys.RESPONSE.getKeyName(), expectedResponse);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    JobService.JobResultOrException<String> resultOrException =
        jobService.retrieveJobResult(flightId, String.class);
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResult_successWithResultTypeRef() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(JobMapKeys.RESPONSE.getKeyName(), expectedResponse);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    JobService.JobResultOrException<String> resultOrException =
        jobService.retrieveJobResult(flightId, null, new TypeReference<>() {});
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResult_fatal() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState fatalFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.FATAL, flightId);
    fatalFlightState.setException(new RuntimeException("test exception"));

    when(mockStairway.getFlightState(any())).thenReturn(fatalFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertEquals(fatalFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResult_error() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState errorFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.ERROR, flightId);
    errorFlightState.setException(new RuntimeException("test exception"));

    when(mockStairway.getFlightState(any())).thenReturn(errorFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertEquals(errorFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResult_running() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState runningFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, flightId);

    when(mockStairway.getFlightState(any())).thenReturn(runningFlightState);

    assertThrows(
        JobNotCompleteException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class, null));
  }

  @Test
  void retrieveJobResult_interrupted() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    when(mockStairway.getFlightState(any())).thenThrow(new InterruptedException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class, null));
  }
}
