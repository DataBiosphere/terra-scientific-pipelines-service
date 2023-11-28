package bio.terra.pipelines.dependencies.common;

public interface HealthCheckWorkspaceApps {

  record Result(boolean isOk, String message) {}

  Result checkHealth(String baseUri, String accessToken);
}
