package bio.terra.pipelines.dependencies.leonardo;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.api.AppsApi;
import org.broadinstitute.dsde.workbench.client.leonardo.api.ServiceInfoApi;
import org.broadinstitute.dsde.workbench.client.leonardo.model.CreateAppRequest;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.broadinstitute.dsde.workbench.client.leonardo.model.SystemStatus;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class LeonardoService implements HealthCheck {

  private final LeonardoClient leonardoClient;
  private final AppUtils appUtils;
  private final RetryTemplate listenerResetRetryTemplate;

  @Autowired
  public LeonardoService(
      LeonardoClient leonardoClient, AppUtils appUtils, RetryTemplate listenerResetRetryTemplate) {
    this.leonardoClient = leonardoClient;
    this.appUtils = appUtils;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  AppsApi getAppsApi(String authToken) {
    return new AppsApi(leonardoClient.getApiClient(authToken));
  }

  ServiceInfoApi getServiceInfoApi() {
    return new ServiceInfoApi(leonardoClient.getUnauthorizedApiClient());
  }

  /** grab app information for a workspace id */
  public List<ListAppResponse> getApps(String workspaceId, String authToken)
      throws LeonardoServiceException {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> getAppsApi(authToken).listAppsV2(workspaceId, null, null, null));
  }

  /**
   * create an app for a given workspace id
   *
   * <p>// // example to create cromwell_runner_app // CreateAppRequest createAppRequest = // new
   * CreateAppRequest() // .appType(AppType.CROMWELL_RUNNER_APP) //
   * .accessScope(AppAccessScope.USER_PRIVATE) // .labels(Map.of("saturnAutoCreated", "true")); //
   * logger.info( // "creating runner app for workspace {}: {}", // workspaceId, //
   * leonardoService.createAppV2( // workspaceId, // samService.getTeaspoonsServiceAccountToken(),
   * // "teasponns-cr-" + UUID.randomUUID(), // createAppRequest));
   *
   * @return
   */
  public Object createAppV2(
      String workspaceId, String authToken, String appName, CreateAppRequest createAppRequest)
      throws LeonardoServiceException {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> {
          getAppsApi(authToken).createAppV2(workspaceId, appName, createAppRequest);
          return null;
        });
  }

  public String getWdsUrlFromApps(String workspaceId, String authToken) {
    return appUtils.findUrlForWds(getApps(workspaceId, authToken), workspaceId);
  }

  public String getCbasUrlFromApps(String workspaceId, String authToken) {
    return appUtils.findUrlForCbas(getApps(workspaceId, authToken), workspaceId);
  }

  public String getCbasUrlFromGetAppResponse(
      List<ListAppResponse> appResponseList, String workspaceId) {
    return appUtils.findUrlForCbas(appResponseList, workspaceId);
  }

  public String getWdsUrlFromGetAppResponse(
      List<ListAppResponse> appResponseList, String workspaceId) {
    return appUtils.findUrlForWds(appResponseList, workspaceId);
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
