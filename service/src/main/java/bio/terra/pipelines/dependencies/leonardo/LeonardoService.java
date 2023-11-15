package bio.terra.pipelines.dependencies.leonardo;

import bio.terra.common.iam.BearerToken;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.api.AppsApi;
import org.broadinstitute.dsde.workbench.client.leonardo.api.ServiceInfoApi;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.broadinstitute.dsde.workbench.client.leonardo.model.SystemStatus;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class LeonardoService implements HealthCheck {

  private final LeonardoClient leonardoClient;
  private final BearerToken bearerToken;
  private final RetryTemplate listenerResetRetryTemplate;

  @Autowired
  public LeonardoService(
      LeonardoClient leonardoClient,
      RetryTemplate listenerResetRetryTemplate,
      BearerToken bearerToken) {
    this.leonardoClient = leonardoClient;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
    this.bearerToken = bearerToken;
  }

  // this will need to be reworked to use the service account credentials instead of the user who
  // made the request
  AppsApi getAppsApi() {
    return new AppsApi(leonardoClient.getApiClient(bearerToken.getToken()));
  }

  ServiceInfoApi getServiceInfoApi() {
    return new ServiceInfoApi(leonardoClient.getUnauthorizedApiClient());
  }

  /** grab app information for a workspace id */
  public List<ListAppResponse> getApps(String workspaceId, boolean creatorOnly)
      throws LeonardoServiceException {
    String creatorRoleSpecifier = creatorOnly ? "creator" : null;
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> getAppsApi().listAppsV2(workspaceId, null, null, null, creatorRoleSpecifier));
  }

  @Override
  public Result checkHealth() {
    try {
      SystemStatus result = getServiceInfoApi().getSystemStatus();
      return new Result(result.getOk(), result.toString());
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  interface LeonardoAction<T> {
    T execute() throws ApiException;
  }

  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, LeonardoAction<T> action)
      throws LeonardoServiceException {

    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new LeonardoServiceApiException(e);
          }
        });
  }
}
