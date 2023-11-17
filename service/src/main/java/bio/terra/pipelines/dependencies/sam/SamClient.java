package bio.terra.pipelines.dependencies.sam;

import bio.terra.common.tracing.OkHttpClientTracingInterceptor;
<<<<<<< HEAD:service/src/main/java/bio/terra/pipelines/dependencies/sam/SamClient.java
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
=======
import bio.terra.pipelines.configuration.external.SamConfiguration;
>>>>>>> 655b57e (folder reconfiguration):service/src/main/java/bio/terra/pipelines/iam/SamClient.java
import io.opencensus.trace.Tracing;
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

  public SamClient(SamConfiguration samConfig) {
    this.samConfig = samConfig;
    this.okHttpClient = new ApiClient().getHttpClient();
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
            .addInterceptor(new OkHttpClientTracingInterceptor(Tracing.getTracer()))
            .build();
    return new ApiClient().setHttpClient(okHttpClientWithTracing).setBasePath(samConfig.basePath());
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
