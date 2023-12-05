package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Optional;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class StairwayJobServiceMockTest extends BaseContainerTest {

  @InjectMocks StairwayJobService stairwayJobService;
  @Mock private Stairway mockStairway;
  @Mock private StairwayComponent mockStairwayComponent;

  @BeforeEach
  void setup() {
    when(mockStairwayComponent.get()).thenReturn(mockStairway);
  }
  /**
   * Reset the {@link StairwayJobService} {@link FlightDebugInfo} after each test so that future
   * submissions aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    stairwayJobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void retrieveJobResult_successWithResultClass() throws InterruptedException {
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(StairwayJobMapKeys.RESPONSE.getKeyName(), expectedResponse);
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatus(FlightStatus.SUCCESS, flightMap);

    String flightId = successFlightState.getFlightId();

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    StairwayJobService.JobResultOrException<String> resultOrException =
        stairwayJobService.retrieveJobResult(flightId, String.class);
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResult_successWithResultTypeRef() throws InterruptedException {
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(StairwayJobMapKeys.RESPONSE.getKeyName(), expectedResponse);
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatus(FlightStatus.SUCCESS, flightMap);

    String flightId = successFlightState.getFlightId();

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    StairwayJobService.JobResultOrException<String> resultOrException =
        stairwayJobService.retrieveJobResult(flightId, null, new TypeReference<>() {});
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResult_fatal() throws InterruptedException {
    FlightState fatalFlightState =
        StairwayTestUtils.constructFlightStateWithStatus(FlightStatus.FATAL);
    fatalFlightState.setException(new RuntimeException("test exception"));
    String flightId = fatalFlightState.getFlightId();

    when(mockStairway.getFlightState(any())).thenReturn(fatalFlightState);

    StairwayJobService.JobResultOrException result =
        stairwayJobService.retrieveJobResult(
            flightId, StairwayJobService.JobResultOrException.class);
    assertEquals(fatalFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResult_running() throws InterruptedException {
    FlightState runningFlightState =
        StairwayTestUtils.constructFlightStateWithStatus(FlightStatus.RUNNING);
    String flightId = runningFlightState.getFlightId();

    when(mockStairway.getFlightState(any())).thenReturn(runningFlightState);

    assertThrows(
        StairwayJobNotCompleteException.class,
        () ->
            stairwayJobService.retrieveJobResult(
                flightId, StairwayJobService.JobResultOrException.class, null));
  }
}
