package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.ForbiddenException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class StairwayJobUnauthorizedException extends ForbiddenException {

  public StairwayJobUnauthorizedException(String message) {
    super(message);
  }

  public StairwayJobUnauthorizedException(String message, Throwable cause) {
    super(message, cause);
  }

  public StairwayJobUnauthorizedException(Throwable cause) {
    super(cause);
  }
}
