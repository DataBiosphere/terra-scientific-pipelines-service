package bio.terra.pipelines.dependencies.cbas;

import bio.terra.common.exception.ErrorReportException;
import java.util.ArrayList;
import org.springframework.http.HttpStatus;

public class CbasServiceException extends ErrorReportException {
  protected CbasServiceException(String message, Throwable cause) {
    super(message, cause, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }

  protected CbasServiceException(String message) {
    super(message, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }
}
