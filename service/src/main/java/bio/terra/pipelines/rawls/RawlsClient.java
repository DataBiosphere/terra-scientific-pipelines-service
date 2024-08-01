package bio.terra.pipelines.rawls;

import bio.terra.catalog.config.RawlsConfiguration;
import bio.terra.common.iam.BearerToken;
import bio.terra.rawls.api.EntitiesApi;
import bio.terra.rawls.api.StatusApi;
import bio.terra.rawls.api.WorkspacesApi;
import bio.terra.rawls.client.ApiClient;
import jakarta.ws.rs.client.Client;
import org.springframework.stereotype.Component;

@Component
public class RawlsClient {
    private final RawlsConfiguration rawlsConfig;
    private final BearerToken bearerToken;
    private final Client commonHttpClient = new ApiClient().getHttpClient();

    public RawlsClient(RawlsConfiguration rawlsConfig, BearerToken bearerToken) {
        this.rawlsConfig = rawlsConfig;
        this.bearerToken = bearerToken;
    }

    private ApiClient getAuthApiClient() {
        ApiClient apiClient = getApiClient();
        apiClient.setAccessToken(bearerToken.getToken());
        return apiClient;
    }

    private ApiClient getApiClient() {
        // Share one api client across requests.
        return new ApiClient().setHttpClient(commonHttpClient).setBasePath(rawlsConfig.basePath());
    }

    WorkspacesApi workspacesApi() {
        return new WorkspacesApi(getAuthApiClient());
    }

    EntitiesApi entitiesApi() {
        return new EntitiesApi(getAuthApiClient());
    }

    StatusApi statusApi() {
        return new StatusApi(getApiClient());
    }
}
