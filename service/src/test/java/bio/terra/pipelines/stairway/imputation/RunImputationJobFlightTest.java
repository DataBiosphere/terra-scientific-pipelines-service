package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class RunImputationJobFlightTest extends BaseEmbeddedDbTest {

  @Autowired private JobService jobService;

  /**
   * How long to wait for a Stairway flight to complete before timing out the test. This is set to 5
   * minutes to allow tests to ride through service outages, cloud retries, and IAM propagation.
   */
  private static final PipelinesEnum imputationPipelineName = PipelinesEnum.IMPUTATION_MINIMAC4;

  private static final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private static final String testUserId = TestUtils.TEST_USER_ID_1;

  private static final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final String testResultPath = TestUtils.TEST_RESULT_URL;

  @Autowired FlightBeanBag flightBeanBag;

  @Test
  void createJobFlightSetup() {
    // this tests the setters for this flight in JobBuilder
    assertDoesNotThrow(
        () ->
            jobService
                .newJob()
                .jobId(testJobId)
                .flightClass(RunImputationJobFlight.class)
                .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), "test RunImputationJobFlight")
                .addParameter(JobMapKeys.USER_ID.getKeyName(), testUserId)
                .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName)
                .addParameter(RunImputationJobFlightMapKeys.PIPELINE_ID, testPipelineId)
                .addParameter(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, testPipelineInputs));
  }

  @Test
  void expectedStepsInFlight() {
    RunImputationJobFlight runImputationJobFlight =
        new RunImputationJobFlight(StairwayTestUtils.CREATE_JOB_INPUT_PARAMS, flightBeanBag);
    assertEquals(8, runImputationJobFlight.getSteps().size());
  }
}
