package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelinesServiceMockTest extends BaseEmbeddedDbTest {
  @Autowired private PipelinesService pipelinesService;
  @MockBean private PipelinesRepository pipelinesRepository;

  @Test
  void getPipelinesOk() {
    // getPipelines should return a list of Pipelines
    List<Pipeline> pipelinesList = List.of(TestUtils.TEST_PIPELINE_1, TestUtils.TEST_PIPELINE_2);
    when(pipelinesRepository.findAll()).thenReturn(pipelinesList);

    List<Pipeline> returnedPipelines = pipelinesService.getPipelines();
    assertEquals(2, returnedPipelines.size());
    assertTrue(returnedPipelines.containsAll(pipelinesList));
  }
}
