package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.testutils.BaseUnitTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class PipelinesServiceTest extends BaseUnitTest {
  @Autowired private PipelinesService pipelinesService;

  @Test
  void testValidatePipeline_exists() {
    // when validating an existing pipeline, should return True
    String existingPipelineId = "imputation";

    // no error should be thrown
    pipelinesService.validatePipeline(existingPipelineId);
  }

  @Test
  void testValidatePipeline_doesNotExist() {
    // when validating an existing pipeline, should return True
    String notExistingPipelineId = "foo";

    assertThrows(
        NotFoundException.class, () -> pipelinesService.validatePipeline(notExistingPipelineId));
  }
}
