package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class RunImputationJobFlightTest extends BaseContainerTest {

  @Autowired private StairwayJobService stairwayJobService;

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final Long STAIRWAY_FLIGHT_TIMEOUT_SECONDS = 300L;

  private static final PipelinesEnum imputationPipelineId = PipelinesEnum.IMPUTATION;
  private static final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private static final String testUserId = TestUtils.TEST_USER_ID_1;

  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;

  @Test
  void createJobFlight_success() throws Exception {
    FlightMap inputParameters =
        StairwayTestUtils.constructCreateJobInputs(
            imputationPipelineId, testPipelineVersion, testUserId, testPipelineInputs);

    FlightState flightState =
        StairwayTestUtils.blockUntilFlightCompletes(
            stairwayJobService.getStairway(),
            RunImputationJobFlight.class,
            UUID.randomUUID(),
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
                .jobId(UUID.randomUUID())
                .flightClass(RunImputationJobFlight.class)
                .pipelineId(imputationPipelineId)
                .pipelineVersion(testPipelineVersion)
                .userId(testUserId)
                .pipelineInputs(testPipelineInputs));
  }
}