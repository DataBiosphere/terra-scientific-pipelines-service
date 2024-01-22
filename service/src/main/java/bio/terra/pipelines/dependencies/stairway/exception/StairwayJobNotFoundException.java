package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.NotFoundException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class StairwayJobNotFoundException extends NotFoundException {
  public StairwayJobNotFoundException(String message) {
    super(message);
  }

  public StairwayJobNotFoundException(String message, Throwable cause) {
    super(message, cause);
  }

  public StairwayJobNotFoundException(Throwable cause) {
    super(cause);
  }
}
