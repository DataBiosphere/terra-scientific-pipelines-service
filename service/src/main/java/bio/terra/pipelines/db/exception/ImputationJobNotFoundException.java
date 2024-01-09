package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.NotFoundException;

public class ImputationJobNotFoundException extends NotFoundException {

  public ImputationJobNotFoundException(String message) {
    super(message);
  }
}
