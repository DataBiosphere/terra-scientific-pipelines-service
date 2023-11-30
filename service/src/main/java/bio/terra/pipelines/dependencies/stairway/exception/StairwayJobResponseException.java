package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.InternalServerErrorException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class StairwayJobResponseException extends InternalServerErrorException {

  public StairwayJobResponseException(String message) {
    super(message);
  }

  public StairwayJobResponseException(String message, Throwable cause) {
    super(message, cause);
  }

  public StairwayJobResponseException(Throwable cause) {
    super(cause);
  }
}
