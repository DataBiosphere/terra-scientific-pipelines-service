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
    PipelineApiUtils.validatePipelineName(PipelinesEnum.IMPUTATION_BEAGLE.getValue(), logger);
  }

  @Test
  void validatePipelineNameBadPipelineName() {
    // give a valid enum and expect no errors
    assertThrows(
        InvalidPipelineException.class,
        () -> PipelineApiUtils.validatePipelineName("bad pipeline name", logger));
  }
}
