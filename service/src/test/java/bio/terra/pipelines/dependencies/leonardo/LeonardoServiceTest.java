package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.net.SocketTimeoutException;
import java.util.*;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.api.AppsApi;
import org.broadinstitute.dsde.workbench.client.leonardo.api.ServiceInfoApi;
import org.broadinstitute.dsde.workbench.client.leonardo.model.CreateAppRequest;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.broadinstitute.dsde.workbench.client.leonardo.model.SystemStatus;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.stubbing.Answer;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

class LeonardoServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks LeonardoService leonardoService;
  @MockBean LeonardoClient leonardoClient;

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

    leonardoService = spy(new LeonardoService(leonardoClient, appUtils, template));
  }

  @Test
  void socketExceptionRetriesEventuallySucceed() throws Exception {
    List<ListAppResponse> expectedResponse = List.of(new ListAppResponse());

    AppsApi appsApi = mock(AppsApi.class);
    when(appsApi.listAppsV2(workspaceId, null, null, null))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(appsApi).when(leonardoService).getAppsApi(authToken);

    assertEquals(expectedResponse, leonardoService.getApps(workspaceId, authToken));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    List<ListAppResponse> expectedResponse = List.of(new ListAppResponse());

    AppsApi appsApi = mock(AppsApi.class);
    when(appsApi.listAppsV2(workspaceId, null, null, null))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(appsApi).when(leonardoService).getAppsApi(authToken);

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          leonardoService.getApps(workspaceId, authToken);
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    List<ListAppResponse> expectedResponse = List.of(new ListAppResponse());

    ApiException expectedException = new ApiException(400, "Bad Leonardo");

    AppsApi appsApi = mock(AppsApi.class);
    when(appsApi.listAppsV2(workspaceId, null, null, null))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    doReturn(appsApi).when(leonardoService).getAppsApi(authToken);

    LeonardoServiceApiException thrown =
        assertThrows(
            LeonardoServiceApiException.class,
            () -> {
              leonardoService.getApps(workspaceId, authToken);
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void getWdsUrlFromApp() {

    List<ListAppResponse> listAppResponseList = List.of(new ListAppResponse().appName("hhaha"));
    doReturn(listAppResponseList).when(leonardoService).getApps(workspaceId, authToken);
    doReturn("bestWdsUri").when(appUtils).findUrlForWds(listAppResponseList, workspaceId);

    assertEquals("bestWdsUri", leonardoService.getWdsUrlFromApps(workspaceId, authToken));
  }

  @Test
  void getWdsUrlFromAppResponse() {
    List<ListAppResponse> listAppResponseList = List.of(new ListAppResponse().appName("hhaasdha"));
    doReturn("bestWdsUri").when(appUtils).findUrlForWds(listAppResponseList, workspaceId);

    assertEquals(
        "bestWdsUri",
        leonardoService.getWdsUrlFromGetAppResponse(listAppResponseList, workspaceId));
  }

  @Test
  void getCbasUrlFromApp() {

    List<ListAppResponse> listAppResponseList = List.of(new ListAppResponse().appName("hhaasdha"));
    doReturn(listAppResponseList).when(leonardoService).getApps(workspaceId, authToken);
    doReturn("bestCbasUri").when(appUtils).findUrlForCbas(listAppResponseList, workspaceId);

    assertEquals("bestCbasUri", leonardoService.getCbasUrlFromApps(workspaceId, authToken));
  }

  @Test
  void getCbasUrlFromAppResponse() {

    List<ListAppResponse> listAppResponseList = List.of(new ListAppResponse().appName("hhaasdha"));
    doReturn("bestCbasUri").when(appUtils).findUrlForCbas(listAppResponseList, workspaceId);

    assertEquals(
        "bestCbasUri",
        leonardoService.getCbasUrlFromGetAppResponse(listAppResponseList, workspaceId));
  }

  @Test
  void createAppV2() {
    Map map = Map.of("blah", "ok");
    String appName = UUID.randomUUID().toString();

    doReturn(map)
        .when(leonardoService)
        .createAppV2(workspaceId, authToken, appName, new CreateAppRequest());

    assertEquals(
        map, leonardoService.createAppV2(workspaceId, authToken, appName, new CreateAppRequest()));
  }

  @Test
  void checkHealth() throws ApiException {
    ServiceInfoApi serviceInfoApi = mock(ServiceInfoApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    when(serviceInfoApi.getSystemStatus()).thenReturn(systemStatus);

    doReturn(serviceInfoApi).when(leonardoService).getServiceInfoApi();
    HealthCheck.Result actualResult = leonardoService.checkHealth();

    assertEquals(
        new HealthCheck.Result(systemStatus.getOk(), systemStatus.toString()), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    ServiceInfoApi serviceInfoApi = mock(ServiceInfoApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(serviceInfoApi.getSystemStatus()).thenThrow(apiException);

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());

    doReturn(serviceInfoApi).when(leonardoService).getServiceInfoApi();
    HealthCheck.Result actualResult = leonardoService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);
  }
}
