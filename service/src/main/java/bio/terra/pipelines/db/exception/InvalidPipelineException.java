package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.BadRequestException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InvalidPipelineException extends BadRequestException {

  public InvalidPipelineException(String message) {
    super(message);
  }
}
