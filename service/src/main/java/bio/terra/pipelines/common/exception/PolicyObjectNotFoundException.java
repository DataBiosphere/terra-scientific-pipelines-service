package bio.terra.pipelines.common.exception;

import bio.terra.common.exception.NotFoundException;

public class PolicyObjectNotFoundException extends NotFoundException {
  public PolicyObjectNotFoundException(String message) {
    super(message);
  }

  public PolicyObjectNotFoundException(String message, Throwable cause) {
    super(message, cause);
  }
}
