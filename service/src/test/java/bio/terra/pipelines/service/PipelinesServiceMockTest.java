package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelinesServiceMockTest extends BaseContainerTest {
  @Autowired private PipelinesService pipelinesService;
  @MockBean private PipelinesRepository pipelinesRepository;

  @Test
  void testGetPipelines() {
    // getPipelines should return a list of Pipelines
    List<Pipeline> pipelinesList = List.of(TestUtils.TEST_PIPELINE_1, TestUtils.TEST_PIPELINE_2);
    when(pipelinesRepository.findAll()).thenReturn(pipelinesList);

    List<Pipeline> returnedPipelines = pipelinesService.getPipelines();
    assertEquals(2, returnedPipelines.size());
    assertTrue(returnedPipelines.containsAll(pipelinesList));
  }

  @Test
  void testGetPipeline_nullResultFromDb() {
    // when retrieving a pipeline that does not exist, should throw an IllegalArgumentException
    String notExistingPipelineId = "notExistingPipeline";
    when(pipelinesRepository.findByPipelineId(notExistingPipelineId)).thenReturn(null);

    Throwable exception =
        assertThrows(
            IllegalArgumentException.class,
            () -> pipelinesService.getPipeline(notExistingPipelineId));
    assertEquals("Pipeline not found for pipelineId notExistingPipeline", exception.getMessage());
  }
}
