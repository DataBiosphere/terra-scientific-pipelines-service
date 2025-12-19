package bio.terra.pipelines.service.exception;

import bio.terra.common.exception.InternalServerErrorException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class PipelineInternalServerException extends InternalServerErrorException {
  public PipelineInternalServerException() {
    super("Something went wrong while running the job. Please reach out to support for help.");
  }
}
