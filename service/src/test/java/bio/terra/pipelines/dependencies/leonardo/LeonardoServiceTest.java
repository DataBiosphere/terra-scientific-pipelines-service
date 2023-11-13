package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.common.iam.BearerToken;
import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import java.net.SocketTimeoutException;
import java.util.List;
import java.util.UUID;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.api.AppsApi;
import org.broadinstitute.dsde.workbench.client.leonardo.api.ServiceInfoApi;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.broadinstitute.dsde.workbench.client.leonardo.model.SystemStatus;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.junit.jupiter.MockitoExtension;
import org.mockito.stubbing.Answer;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

@ExtendWith(MockitoExtension.class)
class LeonardoServiceTest {
  final String workspaceId = UUID.randomUUID().toString();
  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  final BearerToken bearerToken = new BearerToken("");

  final Answer<Object> errorAnswer =
      invocation -> {
        throw new SocketTimeoutException("Timeout");
      };

  @BeforeEach
  void init() {
    FixedBackOffPolicy smallerBackoff = new FixedBackOffPolicy();
    smallerBackoff.setBackOffPeriod(5L); // 5 ms
    template.setBackOffPolicy(smallerBackoff);
  }

  @Test
  void socketExceptionRetriesEventuallySucceed() throws Exception {
    List<ListAppResponse> expectedResponse = List.of(new ListAppResponse());

    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    AppsApi appsApi = mock(AppsApi.class);
    when(appsApi.listAppsV2(workspaceId, null, null, null, null))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    LeonardoService leonardoService =
        spy(new LeonardoService(leonardoClient, template, bearerToken));

    doReturn(appsApi).when(leonardoService).getAppsApi();

    assertEquals(expectedResponse, leonardoService.getApps(workspaceId, false));
  }

  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    List<ListAppResponse> expectedResponse = List.of(new ListAppResponse());

    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    AppsApi appsApi = mock(AppsApi.class);
    when(appsApi.listAppsV2(workspaceId, null, null, null, null))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    LeonardoService leonardoService =
        spy(new LeonardoService(leonardoClient, template, bearerToken));

    doReturn(appsApi).when(leonardoService).getAppsApi();

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          leonardoService.getApps(workspaceId, false);
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    List<ListAppResponse> expectedResponse = List.of(new ListAppResponse());

    ApiException expectedException = new ApiException(400, "Bad Leonardo");

    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    AppsApi appsApi = mock(AppsApi.class);
    when(appsApi.listAppsV2(workspaceId, null, null, null, null))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    LeonardoService leonardoService =
        spy(new LeonardoService(leonardoClient, template, bearerToken));

    doReturn(appsApi).when(leonardoService).getAppsApi();

    LeonardoServiceApiException thrown =
        assertThrows(
            LeonardoServiceApiException.class,
            () -> {
              leonardoService.getApps(workspaceId, false);
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void checkHealth() throws ApiException {
    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    ServiceInfoApi serviceInfoApi = mock(ServiceInfoApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    when(serviceInfoApi.getSystemStatus()).thenReturn(systemStatus);

    LeonardoService leonardoService =
        spy(new LeonardoService(leonardoClient, template, bearerToken));

    doReturn(serviceInfoApi).when(leonardoService).getServiceInfoApi();
    HealthCheck.Result actualResult = leonardoService.checkHealth();

    assertEquals(
        new HealthCheck.Result(systemStatus.getOk(), systemStatus.toString()), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    ServiceInfoApi serviceInfoApi = mock(ServiceInfoApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(serviceInfoApi.getSystemStatus()).thenThrow(apiException);

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());
    LeonardoService leonardoService =
        spy(new LeonardoService(leonardoClient, template, bearerToken));

    doReturn(serviceInfoApi).when(leonardoService).getServiceInfoApi();
    HealthCheck.Result actualResult = leonardoService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);
  }
}
