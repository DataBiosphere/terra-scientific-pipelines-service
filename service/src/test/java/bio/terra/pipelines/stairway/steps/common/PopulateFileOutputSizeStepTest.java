package bio.terra.pipelines.stairway.steps.common;

import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class PopulateFileOutputSizeStepTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesService pipelinesService;
  @Autowired private PipelineRunsRepository pipelineRunsRepository;
  @Autowired private PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired private GcsService gcsService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private FlightMap inputParameters;
  private FlightMap workingMap;

  @BeforeEach
  void setup() {
    inputParameters = new FlightMap();
    workingMap = new FlightMap();

    inputParameters.put(JobMapKeys.PIPELINE_ID, testJobId);
    workingMap.put(
        ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    StairwayTestUtils.constructCreateJobInputs(flightContext.getInputParameters());

    // Todo - add remaining code ....

  }
}
