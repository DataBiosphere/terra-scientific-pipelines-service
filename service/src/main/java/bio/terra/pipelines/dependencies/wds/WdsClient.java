package bio.terra.pipelines.dependencies.wds;

import bio.terra.pipelines.app.configuration.external.WdsServerConfiguration;
import bio.terra.pipelines.dependencies.common.DependencyNotAvailableException;
import okhttp3.OkHttpClient;
import org.databiosphere.workspacedata.api.GeneralWdsInformationApi;
import org.databiosphere.workspacedata.api.RecordsApi;
import org.databiosphere.workspacedata.api.SchemaApi;
import org.databiosphere.workspacedata.client.ApiClient;
import org.springframework.stereotype.Component;

@Component
public class WdsClient {

  private final OkHttpClient singletonHttpClient;
  private final WdsServerConfiguration wdsServerConfiguration;

  public WdsClient(WdsServerConfiguration wdsServerConfiguration) {
    this.wdsServerConfiguration = wdsServerConfiguration;
    singletonHttpClient = new ApiClient().getHttpClient().newBuilder().build();
  }

  protected ApiClient getApiClient(String wdsBaseUri, String accessToken)
      throws DependencyNotAvailableException {

    ApiClient apiClient = new ApiClient().setBasePath(wdsBaseUri);
    apiClient.setHttpClient(singletonHttpClient);
    apiClient.addDefaultHeader("Authorization", "Bearer " + accessToken);
    // By closing the connection after each request, we avoid the problem of the open connection
    // being force-closed ungracefully by the Azure Relay/Listener infrastructure:
    apiClient.addDefaultHeader("Connection", "close");
    apiClient.setDebugging(wdsServerConfiguration.debugApiLogging());
    return apiClient;
  }

  RecordsApi recordsApi(String wdsBaseUri, String accessToken)
      throws DependencyNotAvailableException {
    return new RecordsApi(getApiClient(wdsBaseUri, accessToken));
  }

  GeneralWdsInformationApi generalWdsInformationApi(String wdsBaseUri, String accessToken)
      throws DependencyNotAvailableException {
    return new GeneralWdsInformationApi(getApiClient(wdsBaseUri, accessToken));
  }

  SchemaApi schemaApi(String wdsBaseUri, String accessToken)
      throws DependencyNotAvailableException {
    return new SchemaApi(getApiClient(wdsBaseUri, accessToken));
  }
}
