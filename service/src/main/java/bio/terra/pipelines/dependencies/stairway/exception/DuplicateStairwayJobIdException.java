package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.ConflictException;

/** An exception indicating a jobId is already in use. Error code is 409 CONFLICT. */
public class DuplicateStairwayJobIdException extends ConflictException {
  public DuplicateStairwayJobIdException(String message, Throwable cause) {
    super(message, cause);
  }
}
