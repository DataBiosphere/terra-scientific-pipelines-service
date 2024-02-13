package bio.terra.pipelines.dependencies.cbas;

import bio.terra.cbas.client.ApiException;
import bio.terra.cbas.model.*;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import java.util.UUID;
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

  /**
   * create a cbas wdl method
   *
   * <p>// //example of creating method using cbas service // PostMethodRequest postMethodRequest =
   * // new PostMethodRequest() // .methodName(pipeline.getWdlMethodName()) //
   * .methodDescription("method description") //
   * .methodSource(PostMethodRequest.MethodSourceEnum.GITHUB) // .methodUrl(pipeline.getWdlUrl()) //
   * .methodVersion("1.0"); // logger.info( // "this is creating a new method in cbas: {}", //
   * cbasService.createMethod( // cbasUri, samService.getTspsServiceAccountToken(),
   * postMethodRequest));
   *
   * @param cbasBaseUri - base uri for cbas
   * @param accessToken - tsps SA access token
   * @param postMethodRequest - request capturing method to be created
   * @return - response containing details of method created
   */
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

  public RunLogResponse getRunsForRunSet(String cbasBaseUri, String accessToken, UUID runSetId) {

    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () -> cbasClient.runsApi(cbasBaseUri, accessToken).getRuns(runSetId));
  }

  @Override
  public Result checkHealth(String cbasBaseUri, String accessToken) {
    try {
      SystemStatus result = cbasClient.publicApi(cbasBaseUri, accessToken).getStatus();
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
