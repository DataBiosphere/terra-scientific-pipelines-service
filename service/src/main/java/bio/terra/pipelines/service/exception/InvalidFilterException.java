package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.BadRequestException;

public class InvalidFilterException extends BadRequestException {
  public InvalidFilterException(String message) {
    super(message);
  }
}
