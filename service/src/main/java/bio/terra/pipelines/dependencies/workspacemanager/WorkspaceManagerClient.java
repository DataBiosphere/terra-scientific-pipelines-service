package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.pipelines.app.configuration.external.WorkspaceManagerServerConfiguration;
import bio.terra.workspace.api.ControlledAzureResourceApi;
import bio.terra.workspace.api.ResourceApi;
import bio.terra.workspace.api.UnauthenticatedApi;
import bio.terra.workspace.client.ApiClient;
import org.springframework.stereotype.Component;

@Component
public class WorkspaceManagerClient {

  private final WorkspaceManagerServerConfiguration WorkspaceManagerServerConfiguration;

  public WorkspaceManagerClient(
      WorkspaceManagerServerConfiguration WorkspaceManagerServerConfiguration) {
    this.WorkspaceManagerServerConfiguration = WorkspaceManagerServerConfiguration;
  }

  public ApiClient getUnauthorizedApiClient() {
    var apiClient = new ApiClient().setBasePath(WorkspaceManagerServerConfiguration.baseUri());
    apiClient.setDebugging(WorkspaceManagerServerConfiguration.debugApiLogging());
    return apiClient;
  }

  protected ApiClient getApiClient(String accessToken) {
    ApiClient apiClient = getUnauthorizedApiClient();
    apiClient.setAccessToken(accessToken);
    return apiClient;
  }

  public UnauthenticatedApi getUnauthenticatedApi() {
    return new UnauthenticatedApi(getUnauthorizedApiClient());
  }

  public ResourceApi getResourceApi(String accessToken) {
    return new ResourceApi(getApiClient(accessToken));
  }

  public ControlledAzureResourceApi getControlledAzureResourceApi(String accessToken) {
    return new ControlledAzureResourceApi(getApiClient(accessToken));
  }
}
