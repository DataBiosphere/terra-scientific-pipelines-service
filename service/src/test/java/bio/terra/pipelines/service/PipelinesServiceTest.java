package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.testutils.BaseUnitTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class PipelinesServiceTest extends BaseUnitTest {
  @Autowired private PipelinesService pipelinesService;

  @Test
  void testValidatePipeline_exists() {
    // when validating an existing pipeline, should return True
    String existingPipelineId = "imputation";

    boolean pipelineExists = pipelinesService.validatePipeline(existingPipelineId);

    assertTrue(pipelineExists);
  }

  @Test
  void testValidatePipeline_doesNotExist() {
    // when validating an existing pipeline, should return True
    String notExistingPipelineId = "foo";

    boolean pipelineExists = pipelinesService.validatePipeline(notExistingPipelineId);

    assertFalse(pipelineExists);
  }
}
