package bio.terra.pipelines.dependencies.leonardo;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiClient;
import org.springframework.stereotype.Component;

@Component
public class LeonardoClient {

  private final LeonardoServerConfiguration leonardoServerConfiguration;

  public LeonardoClient(LeonardoServerConfiguration leonardoServerConfiguration) {
    this.leonardoServerConfiguration = leonardoServerConfiguration;
  }

  public ApiClient getUnauthorizedApiClient() {
    var apiClient = new ApiClient().setBasePath(leonardoServerConfiguration.baseUri());
    apiClient.setDebugging(leonardoServerConfiguration.debugApiLogging());
    return apiClient;
  }

  public ApiClient getApiClient(String accessToken) {
    ApiClient apiClient = getUnauthorizedApiClient();
    apiClient.setAccessToken(accessToken);
    return apiClient;
  }
}
