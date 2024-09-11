package bio.terra.pipelines.dependencies.rawls;

import bio.terra.pipelines.app.configuration.external.RawlsConfiguration;
import bio.terra.rawls.api.*;
import bio.terra.rawls.client.ApiClient;
import jakarta.ws.rs.client.Client;
import org.springframework.stereotype.Component;

@Component
public class RawlsClient {
  private final RawlsConfiguration rawlsConfig;
  private final Client commonHttpClient = new ApiClient().getHttpClient();

  public RawlsClient(RawlsConfiguration rawlsConfig) {
    this.rawlsConfig = rawlsConfig;
  }

  protected ApiClient getApiClient(String accessToken) {
    ApiClient apiClient = getUnauthorizedApiClient();
    apiClient.setAccessToken(accessToken);
    return apiClient;
  }

  protected ApiClient getUnauthorizedApiClient() {
    // Share one api client across requests.
    return new ApiClient()
        .setHttpClient(commonHttpClient)
        .setBasePath(rawlsConfig.baseUri())
        .setDebugging(rawlsConfig.debugApiLogging());
  }

  EntitiesApi getEntitiesApi(String accessToken) {
    return new EntitiesApi(getApiClient(accessToken));
  }

  SubmissionsApi getSubmissionsApi(String accessToken) {
    return new SubmissionsApi(getApiClient(accessToken));
  }

  WorkspacesApi getWorkspacesApi(String accessToken) {
    return new WorkspacesApi(getApiClient(accessToken));
  }

  MethodconfigsApi getMethodConfigsApi(String accessToken) {
    return new MethodconfigsApi(getApiClient(accessToken));
  }

  StatusApi getStatusApi() {
    return new StatusApi(getUnauthorizedApiClient());
  }
}
