package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.NotFoundException;

public class PipelineNotFoundException extends NotFoundException {

  public PipelineNotFoundException(String message) {
    super(message);
  }
}
