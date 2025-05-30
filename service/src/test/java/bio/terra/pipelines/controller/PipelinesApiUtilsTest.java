package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.app.controller.PipelineApiUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

class PipelinesApiUtilsTest {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesApiUtilsTest.class);

  @Test
  void validatePipelineNameHappy() {
    // give a valid enum and expect no errors
    PipelineApiUtils.validatePipelineName(PipelinesEnum.ARRAY_IMPUTATION.getValue(), logger);
  }

  @Test
  void validatePipelineNameCaseInsensitive() {
    // give a valid enum and expect no errors
    PipelineApiUtils.validatePipelineName("ARRAY_imputation", logger);
  }

  @Test
  void validatePipelineNameBadPipelineName() {
    // give a bad enum and expect an error
    assertThrows(
        InvalidPipelineException.class,
        () -> PipelineApiUtils.validatePipelineName("bad pipeline name", logger));
  }
}
