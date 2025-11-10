package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.ErrorReportException;
import java.util.ArrayList;
import org.springframework.http.HttpStatus;

public class InvalidFilterException extends ErrorReportException {
  public InvalidFilterException(String message) {
    super(message, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }
}
