package bio.terra.pipelines.dependencies.rawls;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.*;
import static org.mockito.Mockito.doThrow;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.api.StatusApi;
import bio.terra.rawls.client.ApiException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

class RawlsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks RawlsService rawlsService;
  @MockBean RawlsClient rawlsClient;

  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  @BeforeEach
  void init() {
    FixedBackOffPolicy smallerBackoff = new FixedBackOffPolicy();
    smallerBackoff.setBackOffPeriod(5L); // 5 ms
    template.setBackOffPolicy(smallerBackoff);
  }

  @Test
  void checkHealth() throws ApiException {
    StatusApi statusApi = mock(StatusApi.class);

    when(rawlsClient.getStatusApi()).thenReturn(statusApi);
    doNothing().when(statusApi).systemStatus(); // Rawls' systemStatus() is a void method

    HealthCheck.Result actualResult = rawlsService.checkHealth();

    assertEquals(new HealthCheck.Result(true, "Rawls is ok"), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    StatusApi statusApi = mock(StatusApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);

    doReturn(statusApi).when(rawlsClient).getStatusApi();
    doThrow(apiException).when(statusApi).systemStatus();

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());

    HealthCheck.Result actualResult = rawlsService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);
  }
}
