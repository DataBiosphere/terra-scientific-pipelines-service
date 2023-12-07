package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.InternalServerErrorException;

// Exception for failures serializing and deserializing exceptions
@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class ExceptionSerializerException extends InternalServerErrorException {
  public ExceptionSerializerException(String message) {
    super(message);
  }

  public ExceptionSerializerException(String message, Throwable cause) {
    super(message, cause);
  }

  public ExceptionSerializerException(Throwable cause) {
    super(cause);
  }
}
