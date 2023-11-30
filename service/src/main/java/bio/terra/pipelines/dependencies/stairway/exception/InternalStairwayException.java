package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.InternalServerErrorException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InternalStairwayException extends InternalServerErrorException {
  public InternalStairwayException(String message) {
    super(message);
  }

  public InternalStairwayException(String message, Throwable cause) {
    super(message, cause);
  }

  public InternalStairwayException(Throwable cause) {
    super(cause);
  }
}
