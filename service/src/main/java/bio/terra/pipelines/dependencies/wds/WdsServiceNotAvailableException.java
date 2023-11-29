package bio.terra.pipelines.dependencies.wds;

import bio.terra.pipelines.dependencies.common.DependencyNotAvailableException;

public class WdsServiceNotAvailableException extends WdsServiceException {
  private final DependencyNotAvailableException exception;

  public WdsServiceNotAvailableException(DependencyNotAvailableException exception) {
    this.exception = exception;
  }

  @Override
  public synchronized DependencyNotAvailableException getCause() {
    return exception;
  }
}
