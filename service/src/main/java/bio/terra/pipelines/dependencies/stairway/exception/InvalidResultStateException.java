package bio.terra.pipelines.dependencies.stairway.exception;

import bio.terra.common.exception.InternalServerErrorException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InvalidResultStateException extends InternalServerErrorException {
  public InvalidResultStateException(String message) {
    super(message);
  }

  public InvalidResultStateException(String message, Throwable cause) {
    super(message, cause);
  }

  public InvalidResultStateException(Throwable cause) {
    super(cause);
  }

  public static InvalidResultStateException noResultMap() {
    return new InvalidResultStateException("No result map returned from flight");
  }
}
