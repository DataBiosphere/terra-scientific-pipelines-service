package bio.terra.pipelines.dependencies.sam;

import bio.terra.common.tracing.OkHttpClientTracingInterceptor;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import io.opentelemetry.api.OpenTelemetry;
import okhttp3.OkHttpClient;
import org.broadinstitute.dsde.workbench.client.sam.ApiClient;
import org.broadinstitute.dsde.workbench.client.sam.api.ResourcesApi;
import org.broadinstitute.dsde.workbench.client.sam.api.StatusApi;
import org.broadinstitute.dsde.workbench.client.sam.api.UsersApi;
import org.springframework.stereotype.Component;

@Component
public class SamClient {
  private final SamConfiguration samConfig;
  private final OkHttpClient okHttpClient;
  private final OpenTelemetry openTelemetry;

  public SamClient(SamConfiguration samConfig, OpenTelemetry openTelemetry) {
    this.samConfig = samConfig;
    this.okHttpClient = new ApiClient().getHttpClient();
    this.openTelemetry = openTelemetry;
  }

  private ApiClient getApiClient(String accessToken) {
    ApiClient apiClient = getApiClient();
    apiClient.setAccessToken(accessToken);
    return apiClient;
  }

  private ApiClient getApiClient() {
    var okHttpClientWithTracing =
        this.okHttpClient
            .newBuilder()
            .addInterceptor(new OkHttpClientTracingInterceptor(openTelemetry))
            .build();
    return new ApiClient().setHttpClient(okHttpClientWithTracing).setBasePath(samConfig.baseUri());
  }

  UsersApi usersApi(String accessToken) {
    return new UsersApi(getApiClient(accessToken));
  }

  ResourcesApi resourcesApi(String accessToken) {
    return new ResourcesApi(getApiClient(accessToken));
  }

  StatusApi statusApi() {
    return new StatusApi(getApiClient());
  }
}
