package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.common.exception.ErrorReportException;
import java.util.ArrayList;
import org.springframework.http.HttpStatus;

public abstract class WorkspaceManagerServiceException extends ErrorReportException {
  protected WorkspaceManagerServiceException(String message, Throwable cause) {
    super(message, cause, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }

  protected WorkspaceManagerServiceException(String message) {
    super(message, new ArrayList<>(), HttpStatus.INTERNAL_SERVER_ERROR);
  }
}
