package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import org.slf4j.Logger;

public class PipelineApiUtils {

  private PipelineApiUtils() {}

  /**
   * Validates that the pipelineName is a valid pipelineName and returns the Enum value for the
   * pipelineName
   *
   * <p>Note that in PipelinesServiceTest, we check that all the pipelines in the enum exist in the
   * pipelines table
   *
   * @param pipelineName the pipelineName to validate
   * @param logger logger to log messages to
   * @return the Enum value for the pipelineName
   * @throws InvalidPipelineException if the pipelineName is not valid
   */
  public static PipelinesEnum validatePipelineName(String pipelineName, Logger logger) {
    try {
      return PipelinesEnum.valueOf(pipelineName.toUpperCase());
    } catch (IllegalArgumentException e) {
      logger.error("Unknown pipeline name {}", pipelineName);
      throw new InvalidPipelineException(
          String.format("%s is not a valid pipelineName", pipelineName));
    }
  }
}
