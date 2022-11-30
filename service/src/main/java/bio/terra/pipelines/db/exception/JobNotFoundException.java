package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.NotFoundException;

public class JobNotFoundException extends NotFoundException {

  public JobNotFoundException(String message) {
    super(message);
  }
}
