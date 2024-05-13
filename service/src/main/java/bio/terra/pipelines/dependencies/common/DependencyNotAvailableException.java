package bio.terra.pipelines.dependencies.common;

import bio.terra.common.exception.ErrorReportException;

public class DependencyNotAvailableException extends ErrorReportException {
  public DependencyNotAvailableException(String errorMessage) {
    super(errorMessage);
  }

  public static DependencyNotAvailableException formatDependencyNotAvailableExceptionHelper(
      String dependency, String context) {
    return new DependencyNotAvailableException(
        "Dependency not available: %s. %s".formatted(dependency, context));
  }
}
