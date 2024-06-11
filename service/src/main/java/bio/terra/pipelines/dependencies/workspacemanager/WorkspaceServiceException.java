package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.common.exception.ErrorReportException;
import java.util.ArrayList;
import org.springframework.http.HttpStatus;

public abstract class WorkspaceServiceException extends ErrorReportException {
  protected WorkspaceServiceException(String message, Throwable cause) {
    super(message, cause, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }

  protected WorkspaceServiceException(String message) {
    super(message, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }
}
