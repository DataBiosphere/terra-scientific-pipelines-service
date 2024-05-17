package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import io.micrometer.core.instrument.Counter;
import io.micrometer.core.instrument.Metrics;
import io.micrometer.core.instrument.simple.SimpleMeterRegistry;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class WriteJobToDbStepTest extends BaseEmbeddedDbTest {

  @Autowired private ImputationService imputationService;
  @Autowired private JobsRepository jobsRepository;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  private SimpleMeterRegistry meterRegistry;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);

    meterRegistry = new SimpleMeterRegistry();
    Metrics.globalRegistry.add(meterRegistry);
  }

  @AfterEach
  void tearDown() {
    meterRegistry.clear();
    Metrics.globalRegistry.clear();
  }

  @Test
  void doStepSuccess() throws InterruptedException {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // do the step
    var writeJobStep = new WriteJobToDbStep(imputationService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // make sure the job was written to the db
    Job writtenJob =
        jobsRepository
            .findJobByJobIdAndUserId(
                testJobId, inputParams.get(JobMapKeys.USER_ID.getKeyName(), String.class))
            .orElseThrow();
    assertEquals(TestUtils.TEST_PIPELINE_ID_1, writtenJob.getPipelineId());
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void undoStepSuccess() throws InterruptedException {
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    var writeJobStep = new WriteJobToDbStep(imputationService);
    var result = writeJobStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    Counter counter = meterRegistry.find("tsps.pipeline.failed.count").counter();
    assertNotNull(counter);
    assertEquals(1, counter.count());
  }
}
