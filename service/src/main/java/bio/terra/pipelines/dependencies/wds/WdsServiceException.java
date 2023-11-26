package bio.terra.pipelines.dependencies.wds;

public abstract class WdsServiceException extends Exception {
  @Override
  public String getMessage() {
    return getCause().getMessage();
  }
}
