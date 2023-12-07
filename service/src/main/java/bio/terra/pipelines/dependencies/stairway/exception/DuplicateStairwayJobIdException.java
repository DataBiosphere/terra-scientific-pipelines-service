package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.ConflictException;

/** An exception indicating a jobId is already in use. Error code is 409 CONFLICT. */
@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class DuplicateStairwayJobIdException extends ConflictException {
  public DuplicateStairwayJobIdException(String message, Throwable cause) {
    super(message, cause);
  }
}
