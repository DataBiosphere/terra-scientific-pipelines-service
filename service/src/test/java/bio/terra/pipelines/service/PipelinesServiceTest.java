package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.db.PipelinesDao;
import bio.terra.pipelines.db.exception.PipelineNotFoundException;
import bio.terra.pipelines.testutils.BaseUnitTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelinesServiceTest extends BaseUnitTest {
  @Autowired private PipelinesService pipelinesService;
  @MockBean private PipelinesDao mockPipelinesDao;

  @Test
  void testValidatePipeline_exists() {
    // when validating an existing pipeline, should not throw an error
    String existingPipelineId = "existingPipeline";
    when(mockPipelinesDao.checkPipelineExists(existingPipelineId)).thenReturn(true);

    // no error should be thrown
    pipelinesService.validatePipeline(existingPipelineId);
  }

  @Test
  void testValidatePipeline_doesNotExist() {
    // when validating a non-existing pipeline, should throw a NotFoundException
    String notExistingPipelineId = "notExistingPipeline";
    when(mockPipelinesDao.checkPipelineExists(notExistingPipelineId)).thenReturn(false);

    assertThrows(
        PipelineNotFoundException.class,
        () -> pipelinesService.validatePipeline(notExistingPipelineId));
  }
}
