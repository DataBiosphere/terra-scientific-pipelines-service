package bio.terra.pipelines.common.exception;

import bio.terra.common.exception.InternalServerErrorException;

public class InternalTpsErrorException extends InternalServerErrorException {

  public InternalTpsErrorException(String message) {
    super(message);
  }
}
