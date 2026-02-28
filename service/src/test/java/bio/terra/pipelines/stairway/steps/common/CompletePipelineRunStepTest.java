package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.testutils.TestUtils.CONTROL_WORKSPACE_ID;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.stream.Collectors;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class CompletePipelineRunStepTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Autowired private PipelineOutputsRepository pipelineOutputsRepository;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final Integer effectiveQuotaConsumed = 500;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    workingMap.put(
        ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE);
    workingMap.put(ImputationJobMapKeys.EFFECTIVE_QUOTA_CONSUMED, effectiveQuotaConsumed);

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
            TestUtils.TEST_USER_1_ID,
            TestUtils.TEST_PIPELINE_ID_1,
            TestUtils.TEST_TOOL_VERSION_1,
            TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
            TestUtils.CONTROL_WORKSPACE_NAME,
            TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
            TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
            null,
            null,
            CommonPipelineRunStatusEnum.SUCCEEDED,
            TestUtils.TEST_PIPELINE_DESCRIPTION_1,
            null,
            null));

    // do the step
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService);
    var result = writeJobStep.doStep(flightContext);

    // get info from the flight context to run checks
    FlightMap inputParams = flightContext.getInputParameters();

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());

    // make sure the run was updated with isSuccess and quotaConsumed
    PipelineRun writtenJob =
        pipelineRunsRepository
            .findByJobIdAndUserId(testJobId, inputParams.get(JobMapKeys.USER_ID, String.class))
            .orElseThrow();
    assertEquals(CommonPipelineRunStatusEnum.SUCCEEDED, writtenJob.getStatus());
    assertEquals(effectiveQuotaConsumed, writtenJob.getQuotaConsumed());
    assertTrue(writtenJob.getStatus().isSuccess());

    List<PipelineOutput> pipelineOutputList =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunsId(writtenJob.getId());
    assertEquals(2, pipelineOutputList.size());

    // create a map for easier assertion
    Map<String, String> outputMap =
        pipelineOutputList.stream()
            .collect(
                Collectors.toMap(PipelineOutput::getOutputName, PipelineOutput::getOutputValue));

    // assert expected outputs were written with correct values to database
    assertEquals(
        "gs://fc-secure-%s/testFileOutputValue".formatted(CONTROL_WORKSPACE_ID),
        outputMap.get("testFileOutputKey"));
    assertEquals("testStringOutputValue", outputMap.get("testStringOutputKey"));
  }

  // do we want to test how the step handles a failure in the service call?

  @Test
  void undoStepSuccess() {
    var writeJobStep = new CompletePipelineRunStep(pipelineRunsService);
    var result = writeJobStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
