package bio.terra.pipelines.db.exception;

import bio.terra.common.exception.NotFoundException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class ImputationJobNotFoundException extends NotFoundException {

  public ImputationJobNotFoundException(String message) {
    super(message);
  }
}
