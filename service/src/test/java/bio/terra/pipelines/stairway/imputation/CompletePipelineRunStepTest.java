package bio.terra.pipelines.stairway.imputation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
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

class CompletePipelineRunStepTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    workingMap.put(RunImputationJobFlightMapKeys.RAW_OUTPUTS_MAP, TestUtils.TEST_PIPELINE_OUTPUTS);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // write the run to the db
    pipelineRunsRepository.save(
        new PipelineRun(
            testJobId,
            TestUtils.TEST_USER_ID_1,
            TestUtils.TEST_PIPELINE_ID_1,
            TestUtils.CONTROL_WORKSPACE_ID,
            CommonPipelineRunStatusEnum.SUCCEEDED.toString(),
            TestUtils.TEST_PIPELINE_DESCRIPTION_1,
            TestUtils.TEST_RESULT_URL));

    // do the step
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // make sure the run was updated with isSuccess and outputs
    PipelineRun writtenJob =
        pipelineRunsRepository
            .findByJobIdAndUserId(
                testJobId, inputParams.get(JobMapKeys.USER_ID.getKeyName(), String.class))
            .orElseThrow();
    assertTrue(writtenJob.getIsSuccess());
    assertEquals(TestUtils.TEST_PIPELINE_OUTPUTS_STRING, writtenJob.getOutput());
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void undoStepSuccess() {
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService);
    var result = writeJobStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
