package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class PipelinesServiceMockTest extends BaseEmbeddedDbTest {
  @Autowired private PipelinesService pipelinesService;
  @MockitoBean private PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;

  @Test
  void getPipelinesOk() {
    List<PipelineRuntimeMetadata> pipelineRuntimeMetadata =
        List.of(TestUtils.updateTestPipelineRuntime1WithTestValues());
    when(pipelineRuntimeMetadataRepository.findAll()).thenReturn(pipelineRuntimeMetadata);

    List<Pipeline> returnedPipelines = pipelinesService.getPipelines(false);
    // catalog has 2 non-hidden versions of array_imputation (v1 and v2); both are returned
    assertEquals(2, returnedPipelines.size());
    // pipelines are sorted descending by version, so v2 comes first
    assertEquals(TestUtils.TEST_PIPELINE_1.getName(), returnedPipelines.get(0).getName());
    assertEquals(2, returnedPipelines.get(0).getVersion());
  }
}
