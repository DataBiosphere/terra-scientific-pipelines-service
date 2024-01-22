package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.NotFoundException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class JobNotFoundException extends NotFoundException {
  public JobNotFoundException(String message) {
    super(message);
  }

  public JobNotFoundException(String message, Throwable cause) {
    super(message, cause);
  }

  public JobNotFoundException(Throwable cause) {
    super(cause);
  }
}
