package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.MockMvcUtils;
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
    List<Pipeline> pipelinesList =
        List.of(MockMvcUtils.TEST_PIPELINE_1, MockMvcUtils.TEST_PIPELINE_2);
    when(pipelinesRepository.findAll()).thenReturn(pipelinesList);

    List<Pipeline> returnedPipelines = pipelinesService.getPipelines();
    assertEquals(2, returnedPipelines.size());
    assertTrue(returnedPipelines.containsAll(pipelinesList));
  }

  @Test
  void testPipelineExists_true() {
    // when validating an existing pipeline, should return true
    String existingPipelineId = MockMvcUtils.TEST_PIPELINE_ID_1;
    when(pipelinesRepository.existsByPipelineId(existingPipelineId)).thenReturn(true);

    assertTrue(pipelinesService.pipelineExists(existingPipelineId));
  }

  @Test
  void testPipelineExists_false() {
    // when validating a non-existing pipeline, should return false
    String notExistingPipelineId = "notExistingPipeline";
    when(pipelinesRepository.existsByPipelineId(notExistingPipelineId)).thenReturn(false);

    assertFalse(pipelinesService.pipelineExists(notExistingPipelineId));
  }
}
