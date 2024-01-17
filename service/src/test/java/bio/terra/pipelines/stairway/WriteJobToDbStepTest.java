package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.CommonJobStatusEnum;
import bio.terra.pipelines.db.entities.ImputationJob;
import bio.terra.pipelines.db.repositories.ImputationJobsRepository;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class WriteJobToDbStepTest extends BaseEmbeddedDbTest {

  @Autowired private ImputationService imputationService;
  @Autowired private ImputationJobsRepository imputationJobsRepository;
  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void writeJob_doStep_success() throws InterruptedException {
    // setup
    String testJobId = UUID.randomUUID().toString();
    when(flightContext.getFlightId()).thenReturn(testJobId);

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());
    flightContext
        .getWorkingMap()
        .put(RunImputationJobFlightMapKeys.STATUS, CommonJobStatusEnum.SUBMITTED.name());

    // do the step
    var writeJobStep = new WriteJobToDbStep(imputationService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // make sure the job was written to the db
    ImputationJob writtenJob =
        imputationJobsRepository
            .findJobByJobIdAndUserId(
                UUID.fromString(testJobId),
                inputParams.get(JobMapKeys.USER_ID.getKeyName(), String.class))
            .orElseThrow();
    assertEquals(TestUtils.TEST_PIPELINE_VERSION_1, writtenJob.getPipelineVersion());
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void writeJob_undoStep_success() throws InterruptedException {
    var writeJobStep = new WriteJobToDbStep(imputationService);
    var result = writeJobStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
