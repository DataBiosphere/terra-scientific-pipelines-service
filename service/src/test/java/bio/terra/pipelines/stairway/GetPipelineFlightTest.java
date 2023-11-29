package bio.terra.pipelines.stairway;

import bio.terra.pipelines.testutils.BaseTest;
import java.time.Duration;

public class GetPipelineFlightTest extends BaseTest {

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final Duration STAIRWAY_FLIGHT_TIMEOUT = Duration.ofMinutes(5);

  //    @Autowired
  //    private WorkspaceService workspaceService;

  // TODO add tests for retries
}
