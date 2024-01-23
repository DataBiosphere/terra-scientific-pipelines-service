package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.InternalServerErrorException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class JobResponseException extends InternalServerErrorException {

  public JobResponseException(String message) {
    super(message);
  }

  public JobResponseException(String message, Throwable cause) {
    super(message, cause);
  }

  public JobResponseException(Throwable cause) {
    super(cause);
  }
}
