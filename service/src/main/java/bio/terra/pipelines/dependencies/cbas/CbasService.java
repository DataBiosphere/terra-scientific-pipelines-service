package bio.terra.pipelines.dependencies.cbas;

import bio.terra.cbas.client.ApiException;
import bio.terra.cbas.model.*;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import java.util.List;
import java.util.Set;
import java.util.UUID;
import java.util.stream.Collectors;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class CbasService implements HealthCheckWorkspaceApps {
  private final CbasClient cbasClient;
  private final RetryTemplate listenerResetRetryTemplate;

  private static final List<RunState> FINAL_RUN_STATES =
      List.of(
          RunState.COMPLETE,
          RunState.CANCELED,
          RunState.PAUSED,
          RunState.EXECUTOR_ERROR,
          RunState.SYSTEM_ERROR,
          RunState.UNKNOWN);

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
   * cbasService.createMethod( // cbasUri, samService.getTeaspoonsServiceAccountToken(),
   * postMethodRequest));
   *
   * @param cbasBaseUri - base uri for cbas
   * @param accessToken - teaspoons SA access token
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

  public static UUID getMethodVersionIdFromMethodListResponse(
      MethodListResponse methodListResponse, String methodName) {
    UUID methodVersionId = null;
    for (MethodDetails methodDetails : methodListResponse.getMethods()) {
      if (methodDetails.getName().equals(methodName)) {
        // for now grabbing the first MethodVersionId but should change once we start having a new
        // pipeline for each version of a wdl.
        methodVersionId = methodDetails.getMethodVersions().get(0).getMethodVersionId();
        break;
      }
    }
    return methodVersionId;
  }

  public static boolean containsRunningRunLog(RunLogResponse runLogResponse) {
    Set<RunState> runningRuns =
        runLogResponse.getRuns().stream()
            .map(RunLog::getState)
            .filter(runState -> !isFinalRunState(runState))
            .collect(Collectors.toSet());
    return !runningRuns.isEmpty();
  }

  private static boolean isFinalRunState(RunState runState) {
    return FINAL_RUN_STATES.contains(runState);
  }
}
