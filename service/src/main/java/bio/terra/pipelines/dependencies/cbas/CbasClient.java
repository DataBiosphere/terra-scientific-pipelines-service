package bio.terra.pipelines.dependencies.cbas;

import bio.terra.cbas.api.MethodsApi;
import bio.terra.cbas.api.PublicApi;
import bio.terra.cbas.api.RunSetsApi;
import bio.terra.cbas.api.RunsApi;
import bio.terra.cbas.client.ApiClient;
import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import org.springframework.stereotype.Component;

@Component
public class CbasClient {

  private final CbasConfiguration cbasConfiguration;

  public CbasClient(CbasConfiguration cbasConfiguration) {
    this.cbasConfiguration = cbasConfiguration;
  }

  protected ApiClient getApiClient(String cbasBaseUri, String accessToken) {
    ApiClient apiClient = new ApiClient().setBasePath(cbasBaseUri);
    apiClient.addDefaultHeader("Authorization", "Bearer " + accessToken);
    // By closing the connection after each request, we avoid the problem of the open connection
    // being force-closed ungracefully by the Azure Relay/Listener infrastructure:
    apiClient.addDefaultHeader("Connection", "close");
    apiClient.setDebugging(cbasConfiguration.getDebugApiLogging());
    return apiClient;
  }

  // used to access public endpoints like status and health
  PublicApi publicApi(String cbasBaseUri, String accessToken) {
    return new PublicApi(getApiClient(cbasBaseUri, accessToken));
  }

  // used to access endpoints related to methods (wdls)
  // create, search, delete,etc.
  MethodsApi methodsApi(String cbasBaseUri, String accessToken) {
    return new MethodsApi(getApiClient(cbasBaseUri, accessToken));
  }

  // used to access endpoints to get information related to runs in a run set
  RunsApi runsApi(String cbasBaseUri, String accessToken) {
    return new RunsApi(getApiClient(cbasBaseUri, accessToken));
  }

  // used to access endpoints related to running of a submission of methods
  // launch, abort, list run sets for a given method
  RunSetsApi runSetsApi(String cbasBaseUri, String accessToken) {
    return new RunSetsApi(getApiClient(cbasBaseUri, accessToken));
  }
}
