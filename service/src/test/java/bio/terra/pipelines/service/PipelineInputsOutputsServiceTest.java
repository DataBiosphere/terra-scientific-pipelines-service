package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.anyString;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputs;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.Entity;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelineInputsOutputsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelineInputsOutputsService pipelineInputsOutputsService;

  @Autowired PipelineInputsRepository pipelineInputsRepository;
  @Autowired PipelineOutputsRepository pipelineOutputsRepository;

  @Autowired PipelineRunsRepository pipelineRunsRepository;

  @MockBean private GcsService mockGcsService;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @Test
  void extractPipelineOutputsFromEntity() {
    // test that the method correctly extracts the outputs from the entity
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(
        Map.of("output_name", "gs://bucket/file1", "testNonOutputKey", "doesn't matter"));

    Map<String, String> extractedOutputs =
        pipelineInputsOutputsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);

    assertEquals(1, extractedOutputs.size());
    // the meethod should also have converted the wdlVariableName key to the camelCase outputName
    // key
    assertEquals("gs://bucket/file1", extractedOutputs.get("outputName"));
  }

  @Test
  void extractPipelineOutputsFromEntityMissingOutput() {
    // test that the method correctly throws an error if an output is missing
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(Map.of("testNonOutputKey", "doesn't matter"));

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineInputsOutputsService.extractPipelineOutputsFromEntity(
                outputDefinitions, entity));
  }

  @Test
  void extractPipelineOutputsFromEntityEmptyOutput() {
    // test that the method correctly throws an error if an output is empty
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(Map.of("outputName", ""));

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineInputsOutputsService.extractPipelineOutputsFromEntity(
                outputDefinitions, entity));
  }

  @Test
  void formatPipelineRunOutputs() throws MalformedURLException {
    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(testJobId);
    pipelineRunsRepository.save(pipelineRun);

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(pipelineRun.getId());
    pipelineOutput.setOutputs(
        pipelineInputsOutputsService.mapToString(TestUtils.TEST_PIPELINE_OUTPUTS));
    pipelineOutputsRepository.save(pipelineOutput);

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");
    // mock GCS service
    when(mockGcsService.generateGetObjectSignedUrl(
            eq(pipelineRun.getWorkspaceGoogleProject()),
            eq(pipelineRun.getWorkspaceStorageContainerName()),
            anyString()))
        .thenReturn(fakeUrl);

    ApiPipelineRunOutputs apiPipelineRunOutputs =
        pipelineInputsOutputsService.formatPipelineRunOutputs(pipelineRun);

    assertEquals(fakeUrl.toString(), apiPipelineRunOutputs.get("testFileOutputKey"));
  }
}
