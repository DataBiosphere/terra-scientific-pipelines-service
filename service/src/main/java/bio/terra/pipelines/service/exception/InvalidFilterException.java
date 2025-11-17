package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.BadRequestException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InvalidFilterException extends BadRequestException {
  public InvalidFilterException(String message) {
    super(message);
  }
}
