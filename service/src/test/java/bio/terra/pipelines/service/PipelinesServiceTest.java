package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.PipelinesDao;
import bio.terra.pipelines.testutils.BaseUnitTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelinesServiceTest extends BaseUnitTest {
  @Autowired private PipelinesService pipelinesService;
  @MockBean private PipelinesDao mockPipelinesDao;

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
