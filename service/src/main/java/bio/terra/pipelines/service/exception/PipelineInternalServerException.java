package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.InternalServerErrorException;

public class PipelineInternalServerException extends InternalServerErrorException {
  public PipelineInternalServerException() {
    super("Something went wrong while running the job. Please reach out to support for help.");
  }
}
