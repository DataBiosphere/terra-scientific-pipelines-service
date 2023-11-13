package bio.terra.pipelines.service.stairwayjob.exception;

import bio.terra.common.exception.InternalServerErrorException;

public class StairwayJobResponseException extends InternalServerErrorException {

  public StairwayJobResponseException(String message) {
    super(message);
  }

  public StairwayJobResponseException(String message, Throwable cause) {
    super(message, cause);
  }

  public StairwayJobResponseException(Throwable cause) {
    super(cause);
  }
}
