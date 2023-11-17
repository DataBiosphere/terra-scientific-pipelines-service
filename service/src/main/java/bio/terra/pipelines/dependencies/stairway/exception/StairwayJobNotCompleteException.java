package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.BadRequestException;

public class StairwayJobNotCompleteException extends BadRequestException {
  public StairwayJobNotCompleteException(String message) {
    super(message);
  }

  public StairwayJobNotCompleteException(String message, Throwable cause) {
    super(message, cause);
  }

  public StairwayJobNotCompleteException(Throwable cause) {
    super(cause);
  }
}
