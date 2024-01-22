package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.ForbiddenException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class JobUnauthorizedException extends ForbiddenException {

  public JobUnauthorizedException(String message) {
    super(message);
  }

  public JobUnauthorizedException(String message, Throwable cause) {
    super(message, cause);
  }

  public JobUnauthorizedException(Throwable cause) {
    super(cause);
  }
}
