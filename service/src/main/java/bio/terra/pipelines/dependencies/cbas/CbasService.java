package bio.terra.pipelines.dependencies.cbas;

import bio.terra.cbas.client.ApiException;
import bio.terra.cbas.model.*;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class CbasService implements HealthCheckWorkspaceApps {
  private final CbasClient cbasClient;
  private final RetryTemplate listenerResetRetryTemplate;

  public CbasService(CbasClient cbasClient, RetryTemplate listenerResetRetryTemplate) {
    this.cbasClient = cbasClient;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  public MethodListResponse getAllMethods(String cbasBaseUri, String accessToken) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> cbasClient.methodsApi(cbasBaseUri, accessToken).getMethods(null, null, null));
  }

  // create a cbas method
  public PostMethodResponse createMethod(
      String cbasBaseUri, String accessToken, PostMethodRequest postMethodRequest) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> cbasClient.methodsApi(cbasBaseUri, accessToken).postMethod(postMethodRequest));
  }

  public RunSetStateResponse createRunSet(
      String cbasBaseUri, String accessToken, RunSetRequest runSetRequest) {

    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> cbasClient.runSetsApi(cbasBaseUri, accessToken).postRunSet(runSetRequest));
  }

  @Override
  public Result checkHealth(String wdsBaseUri, String accessToken) {
    try {
      SystemStatus result = cbasClient.publicApi(wdsBaseUri, accessToken).getStatus();
      return new Result(result.isOk(), result.toString());
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  interface CbasAction<T> {
    T execute() throws ApiException;
  }

  static <T> T executionWithRetryTemplate(
      RetryTemplate retryTemplate, CbasService.CbasAction<T> action)
      throws CbasServiceApiException {

    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new CbasServiceApiException(e);
          }
        });
  }
}
