package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.BadRequestException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class JobNotCompleteException extends BadRequestException {
  public JobNotCompleteException(String message) {
    super(message);
  }

  public JobNotCompleteException(String message, Throwable cause) {
    super(message, cause);
  }

  public JobNotCompleteException(Throwable cause) {
    super(cause);
  }
}
