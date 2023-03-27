package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.PipelinesDao;
import bio.terra.pipelines.service.model.Pipeline;
import bio.terra.pipelines.testutils.BaseUnitTest;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelinesServiceTest extends BaseUnitTest {
  @Autowired private PipelinesService pipelinesService;
  @MockBean private PipelinesDao mockPipelinesDao;

  @Test
  void testGetPipelines() {
    // getPipelines should return a list of Pipelines
    Pipeline testPipeline1 =
        new Pipeline("testPipeline1", "Test Pipeline One", "Test Pipeline One Description");
    Pipeline testPipeline2 =
        new Pipeline("testPipeline2", "Test Pipeline Two", "Test Pipeline Two Description");
    List<Pipeline> pipelinesList = List.of(testPipeline1, testPipeline2);
    when(mockPipelinesDao.getPipelines()).thenReturn(pipelinesList);

    List<Pipeline> returnedPipelines = pipelinesService.getPipelines();
    assertEquals(2, returnedPipelines.size());
    assertTrue(returnedPipelines.containsAll(pipelinesList));
  }

  @Test
  void testPipelineExists_true() {
    // when validating an existing pipeline, should return true
    String existingPipelineId = "existingPipeline";
    when(mockPipelinesDao.checkPipelineExists(existingPipelineId)).thenReturn(true);

    assertTrue(pipelinesService.pipelineExists(existingPipelineId));
  }

  @Test
  void testPipelineExists_false() {
    // when validating a non-existing pipeline, should return false
    String notExistingPipelineId = "notExistingPipeline";
    when(mockPipelinesDao.checkPipelineExists(notExistingPipelineId)).thenReturn(false);

    assertFalse(pipelinesService.pipelineExists(notExistingPipelineId));
  }
}
