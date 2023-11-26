package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import java.net.SocketTimeoutException;
import java.util.Collections;
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

  final AppUtils appUtils = mock(AppUtils.class);

  final String authToken = "authToken";

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

    LeonardoService leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));

    doReturn(appsApi).when(leonardoService).getAppsApi(authToken);

    assertEquals(expectedResponse, leonardoService.getApps(workspaceId, authToken, false));
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

    LeonardoService leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));

    doReturn(appsApi).when(leonardoService).getAppsApi(authToken);

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          leonardoService.getApps(workspaceId, authToken, false);
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

    LeonardoService leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));

    doReturn(appsApi).when(leonardoService).getAppsApi(authToken);

    LeonardoServiceApiException thrown =
        assertThrows(
            LeonardoServiceApiException.class,
            () -> {
              leonardoService.getApps(workspaceId, authToken, false);
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void getWdsUrlFromApp() {
    LeonardoClient leonardoClient = mock(LeonardoClient.class);

    LeonardoService leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));

    doReturn(Collections.emptyList()).when(leonardoService).getApps(workspaceId, authToken, false);
    doReturn("bestWdsUri").when(appUtils).findUrlForWds(any(), any());

    assertEquals("bestWdsUri", leonardoService.getWdsUrlFromApps(workspaceId, authToken, false));
  }

  @Test
  void checkHealth() throws ApiException {
    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    ServiceInfoApi serviceInfoApi = mock(ServiceInfoApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    when(serviceInfoApi.getSystemStatus()).thenReturn(systemStatus);

    LeonardoService leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));

    doReturn(serviceInfoApi).when(leonardoService).getServiceInfoApi();
    HealthCheck.Result actualResult = leonardoService.checkHealth();

    assertEquals(
        new HealthCheck.Result(systemStatus.getOk(), systemStatus.toString()), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    LeonardoClient leonardoClient = mock(LeonardoClient.class);
    ServiceInfoApi serviceInfoApi = mock(ServiceInfoApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(serviceInfoApi.getSystemStatus()).thenThrow(apiException);

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());
    LeonardoService leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));

    doReturn(serviceInfoApi).when(leonardoService).getServiceInfoApi();
    HealthCheck.Result actualResult = leonardoService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);
  }
}
