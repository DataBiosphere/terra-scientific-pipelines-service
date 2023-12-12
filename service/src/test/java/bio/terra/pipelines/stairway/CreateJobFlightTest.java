package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class CreateJobFlightTest extends BaseContainerTest {

  @Autowired private StairwayJobService stairwayJobService;

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final Long STAIRWAY_FLIGHT_TIMEOUT_SECONDS = 300L;

  private static final String testPipelineId = MockMvcUtils.TEST_PIPELINE_ID_1;
  private static final String testPipelineVersion = MockMvcUtils.TEST_PIPELINE_VERSION_1;
  private static final String testUserId = MockMvcUtils.TEST_USER_ID_1;
  private static final String testJobId = MockMvcUtils.TEST_UUID_STRING;

  private final Object testPipelineInputs = MockMvcUtils.TEST_PIPELINE_INPUTS;

  @Test
  void createJobFlight_success() throws Exception {
    FlightMap inputParameters =
        StairwayTestUtils.constructCreateJobInputs(
            testPipelineId, testPipelineVersion, testUserId, testPipelineInputs);

    FlightState flightState =
        StairwayTestUtils.blockUntilFlightCompletes(
            stairwayJobService.getStairway(),
            CreateJobFlight.class,
            testJobId,
            inputParameters,
            STAIRWAY_FLIGHT_TIMEOUT_SECONDS,
            /*debugInfo*/ null);

    assertEquals(FlightStatus.SUCCESS, flightState.getFlightStatus());
  }

  @Test
  void createJobFlight_setup() {
    // this tests the setters for this flight in StairwayJobBuilder
    assertDoesNotThrow(
        () ->
            stairwayJobService
                .newJob()
                .jobId(testJobId)
                .flightClass(CreateJobFlight.class)
                .pipelineId(testPipelineId)
                .pipelineVersion(testPipelineVersion)
                .submittingUserId(testUserId)
                .pipelineInputs(testPipelineInputs));
  }
}
