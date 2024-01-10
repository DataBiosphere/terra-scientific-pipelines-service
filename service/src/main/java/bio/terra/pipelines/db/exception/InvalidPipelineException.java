package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.BadRequestException;

public class InvalidPipelineException extends BadRequestException {

  public InvalidPipelineException(String message) {
    super(message);
  }
}
