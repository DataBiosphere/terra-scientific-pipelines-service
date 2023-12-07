package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.BadRequestException;

/** An exception indicating an invalid jobId string value. Error code is 400 Bad Request. */
@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InvalidStairwayJobIdException extends BadRequestException {
  public InvalidStairwayJobIdException(String message) {
    super(message);
  }

  public InvalidStairwayJobIdException(String message, Throwable cause) {
    super(message, cause);
  }

  public InvalidStairwayJobIdException(Throwable cause) {
    super(cause);
  }
}
