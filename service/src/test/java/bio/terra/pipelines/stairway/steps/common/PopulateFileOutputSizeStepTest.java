package bio.terra.pipelines.stairway.steps.common;

import static bio.terra.pipelines.testutils.TestUtils.*;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.stairway.flights.wdlbasedpipelinerun.WdlBasedPipelineJobMapKeys;
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

  @Autowired private PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired private PipelineConfigurations pipelineConfigurations;
  @MockitoBean private GcsService gcsService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();

    inputParameters.put(JobMapKeys.PIPELINE_KEY, TestUtils.TEST_PIPELINE_KEY_1);
    workingMap.put(
        WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS,
        new HashMap<>(
            Map.of(
                "testOutput",
                "gs://fc-secure-%s/test-output.vcf.gz".formatted(CONTROL_WORKSPACE_ID),
                "testStringOutputKey",
                "not-a-file-output")));

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccessWithFileSizes() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    when(gcsService.getFileSizeInBytes(
            "gs://fc-secure-%s/test-output.vcf.gz".formatted(CONTROL_WORKSPACE_ID)))
        .thenReturn(256L);

    // do the step
    var populateFileSizeStep = new PopulateFileOutputSizeStep(pipelineInputsOutputsService);
    var result = populateFileSizeStep.doStep(flightContext);

    // verify step success and that file sizes were populated in the working map correctly
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        Map.of("testOutput", 256L),
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class));
  }

  @Test
  void doStepSuccessWithoutFileSizes() {
    // setup
    when(flightContext.getFlightId()).thenReturn(testJobId.toString());

    // simulate an error when trying to get the file size from GCS, which should be caught
    // and logged in the step, and the step should continue without populating file sizes
    // in the working map
    when(gcsService.getFileSizeInBytes(
            "gs://fc-secure-%s/test-output.vcf.gz".formatted(CONTROL_WORKSPACE_ID)))
        .thenThrow(new InternalServerErrorException("GCS Service Exception thrown for testing"));

    // do the step
    var populateFileSizeStep = new PopulateFileOutputSizeStep(pipelineInputsOutputsService);
    var result = populateFileSizeStep.doStep(flightContext);

    // verify step success and that file sizes were not populated in the working map
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertNull(
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class));
  }

  @Test
  void undoStepSuccess() {
    // do the undo
    var populateFileSizeStep = new PopulateFileOutputSizeStep(pipelineInputsOutputsService);
    var result = populateFileSizeStep.undoStep(flightContext);

    // verify step success
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }

  @Test
  void doStepUsesPipelineKeyYamlDefinitions() {
    FlightMap inputParameters = new FlightMap();
    FlightMap workingMap = new FlightMap();
    inputParameters.put(JobMapKeys.PIPELINE_KEY, "array_imputation_v1");
    workingMap.put(
        WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS,
        new HashMap<>(
            Map.of(
                "testOutput",
                "gs://fc-secure-%s/test-output.vcf.gz".formatted(CONTROL_WORKSPACE_ID),
                "testStringOutputKey",
                "not-a-file-output")));

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
    when(gcsService.getFileSizeInBytes(
            "gs://fc-secure-%s/test-output.vcf.gz".formatted(CONTROL_WORKSPACE_ID)))
        .thenReturn(256L);

    var populateFileSizeStep = new PopulateFileOutputSizeStep(pipelineInputsOutputsService);
    var result = populateFileSizeStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    assertEquals(
        Map.of("testOutput", 256L),
        flightContext
            .getWorkingMap()
            .get(WdlBasedPipelineJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, Map.class));
  }
}
