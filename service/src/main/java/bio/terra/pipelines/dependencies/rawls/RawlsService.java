package bio.terra.pipelines.dependencies.rawls;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.rawls.client.ApiException;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class RawlsService implements HealthCheck {
  private final RawlsClient rawlsClient;
  private final RetryTemplate listenerResetRetryTemplate;

  public RawlsService(RawlsClient rawlsClient, RetryTemplate listenerResetRetryTemplate) {
    this.rawlsClient = rawlsClient;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  @Override
  public Result checkHealth() {
    // systemStatus is a void method, throws ApiException if status is not healthy
    try {
      rawlsClient.getStatusApi().systemStatus();
      return new Result(true, "Rawls is ok");
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  interface RawlsAction<T> {
    T execute() throws ApiException;
  }

  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, RawlsAction<T> action)
      throws RawlsServiceApiException {

    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new RawlsServiceApiException(e);
          }
        });
  }
}
