package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.PipelinesDao;
import bio.terra.pipelines.service.model.Pipeline;
import bio.terra.pipelines.testutils.BaseUnitTest;
import bio.terra.pipelines.testutils.MockMvcUtils;
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
    List<Pipeline> pipelinesList =
        List.of(MockMvcUtils.TEST_PIPELINE_1, MockMvcUtils.TEST_PIPELINE_2);
    when(mockPipelinesDao.getPipelines()).thenReturn(pipelinesList);

    List<Pipeline> returnedPipelines = pipelinesService.getPipelines();
    assertEquals(2, returnedPipelines.size());
    assertTrue(returnedPipelines.containsAll(pipelinesList));
  }

  @Test
  void testPipelineExists_true() {
    // when validating an existing pipeline, should return true
    String existingPipelineId = MockMvcUtils.TEST_PIPELINE_ID_1;
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
