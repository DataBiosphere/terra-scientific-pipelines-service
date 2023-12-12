package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import java.util.HashMap;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class CreateJobFlightTest extends BaseContainerTest {

  @Autowired private StairwayJobService stairwayJobService;

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final Long STAIRWAY_FLIGHT_TIMEOUT_SECONDS = 300L;

  @Test
  void createJobFlight_success() throws Exception {
    String pipelineId = "imputation";
    String pipelineVersion = "v0";
    String submittingUserId = "submittingUserId";
    Object pipelineInputs = new HashMap<>();
    FlightMap inputParameters =
        StairwayTestUtils.constructCreateJobInputs(
            pipelineId, pipelineVersion, submittingUserId, pipelineInputs);

    FlightState flightState =
        StairwayTestUtils.blockUntilFlightCompletes(
            stairwayJobService.getStairway(),
            CreateJobFlight.class,
            UUID.randomUUID().toString(),
            inputParameters,
            STAIRWAY_FLIGHT_TIMEOUT_SECONDS,
            /*debugInfo*/ null);

    assertEquals(FlightStatus.SUCCESS, flightState.getFlightStatus());
  }
}
