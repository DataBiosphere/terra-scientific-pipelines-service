package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.testutils.TestUtils.*;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class PopulateFileOutputSizeStepTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesService pipelinesService;
  @Autowired private PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired private PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;
  @MockitoBean private GcsService gcsService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private Map<String, String> outputValuesByName;
  private Map<String, Long> expectedOutputFileSizes;

  @BeforeEach
  void setup() {
    PipelineRuntimeMetadata testPipeline =
        pipelineRuntimeMetadataRepository.findByNameAndVersion(PipelinesEnum.ARRAY_IMPUTATION, 1);
    var configuredPipeline = pipelinesService.getPipelineById(testPipeline.getId());

    Long testPipelineId = testPipeline.getId();

    // set up the flight context input and working maps
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    outputValuesByName = new HashMap<>();
    expectedOutputFileSizes = new HashMap<>();
    for (PipelineOutputDefinition outputDefinition :
        configuredPipeline.getPipelineOutputDefinitions()) {
      if (outputDefinition.getType() == PipelineVariableTypesEnum.FILE) {
        String outputPath =
            "gs://fc-secure-%s/%s.vcf.gz"
                .formatted(CONTROL_WORKSPACE_ID, outputDefinition.getName());
        outputValuesByName.put(outputDefinition.getName(), outputPath);
        expectedOutputFileSizes.put(outputDefinition.getName(), 256L);
      }
    }

    inputParameters.put(JobMapKeys.PIPELINE_ID, testPipelineId);
    workingMap.put(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, outputValuesByName);

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccessWithFileSizes() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    outputValuesByName
        .values()
        .forEach(path -> when(gcsService.getFileSizeInBytes(path)).thenReturn(256L));

    // do the step
    var populateFileSizeStep =
        new PopulateFileOutputSizeStep(pipelinesService, pipelineInputsOutputsService);
    var result = populateFileSizeStep.doStep(flightContext);

    // verify step success and that file sizes were populated in the working map correctly
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        expectedOutputFileSizes,
        flightContext
            .getWorkingMap()
            .get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class));
  }

  @Test
  void doStepSuccessWithoutFileSizes() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    // simulate an error when trying to get the file size from GCS, which should be caught
    // and logged in the step, and the step should continue without populating file sizes
    // in the working map
    when(gcsService.getFileSizeInBytes(outputValuesByName.values().iterator().next()))
        .thenThrow(new InternalServerErrorException("GCS Service Exception thrown for testing"));

    // do the step
    var populateFileSizeStep =
        new PopulateFileOutputSizeStep(pipelinesService, pipelineInputsOutputsService);
    var result = populateFileSizeStep.doStep(flightContext);

    // verify step success and that file sizes were not populated in the working map
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertNull(
        flightContext
            .getWorkingMap()
            .get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class));
  }

  @Test
  void undoStepSuccess() {
    // do the undo
    var populateFileSizeStep =
        new PopulateFileOutputSizeStep(pipelinesService, pipelineInputsOutputsService);
    var result = populateFileSizeStep.undoStep(flightContext);

    // verify step success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
