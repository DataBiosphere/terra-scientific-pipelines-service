package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.ErrorReportException;

public class PipelineCheckFailedException extends ErrorReportException {
  public PipelineCheckFailedException(String message) {
    super(message);
  }
}
