package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.BadRequestException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class RequesterPaysBucketException extends BadRequestException {
  public RequesterPaysBucketException(String message) {
    super(message);
  }
}
