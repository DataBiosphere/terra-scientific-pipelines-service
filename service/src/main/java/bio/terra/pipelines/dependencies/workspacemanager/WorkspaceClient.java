package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.pipelines.app.configuration.external.WorkspaceServerConfiguration;
import bio.terra.workspace.client.ApiClient;
import org.springframework.stereotype.Component;

@Component
public class WorkspaceClient {

  private final WorkspaceServerConfiguration workspaceServerConfiguration;

  public WorkspaceClient(WorkspaceServerConfiguration workspaceServerConfiguration) {
    this.workspaceServerConfiguration = workspaceServerConfiguration;
  }

  public ApiClient getUnauthorizedApiClient() {
    var apiClient = new ApiClient().setBasePath(workspaceServerConfiguration.baseUri());
    apiClient.setDebugging(workspaceServerConfiguration.debugApiLogging());
    return apiClient;
  }

  protected ApiClient getApiClient(String accessToken) {
    ApiClient apiClient = getUnauthorizedApiClient();
    apiClient.setAccessToken(accessToken);
    return apiClient;
  }
}
