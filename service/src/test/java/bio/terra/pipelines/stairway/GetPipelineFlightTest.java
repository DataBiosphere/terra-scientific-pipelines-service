package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import java.time.Duration;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class GetPipelineFlightTest extends BaseContainerTest {

  @Autowired private StairwayJobService stairwayJobService;

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final Duration STAIRWAY_FLIGHT_TIMEOUT = Duration.ofMinutes(5);

  @Test
  void getPipelineFlight_success() throws Exception {
    String pipelineId = "imputation";
    FlightMap inputParameters = StairwayTestUtils.constructPipelineInputs(pipelineId);

    FlightState flightState =
        StairwayTestUtils.blockUntilFlightCompletes(
            stairwayJobService.getStairway(),
            GetPipelineFlight.class,
            inputParameters,
            STAIRWAY_FLIGHT_TIMEOUT,
            /*debugInfo*/ null);

    assertEquals(FlightStatus.SUCCESS, flightState.getFlightStatus());
  }
}
