package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.BadRequestException;

public class DuplicateObjectException extends BadRequestException {
  public DuplicateObjectException(String message) {
    super(message);
  }

  public DuplicateObjectException(String message, Throwable cause) {
    super(message, cause);
  }
}
