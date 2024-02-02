package bio.terra.pipelines.dependencies.cbas;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;
import static org.mockito.Mockito.doReturn;

import bio.terra.cbas.api.MethodsApi;
import bio.terra.cbas.api.PublicApi;
import bio.terra.cbas.api.RunSetsApi;
import bio.terra.cbas.client.ApiException;
import bio.terra.cbas.model.*;
import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import java.net.SocketTimeoutException;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.junit.jupiter.MockitoExtension;
import org.mockito.stubbing.Answer;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

@ExtendWith(MockitoExtension.class)
class CbasServiceTest {
  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  final CbasConfiguration cbasConfiguration = new CbasConfiguration(false);

  final Answer<Object> errorAnswer =
      invocation -> {
        throw new SocketTimeoutException("Timeout");
      };

  private final String cbaseBaseUri = "baseUri";
  private final String accessToken = "accessToken";

  @BeforeEach
  void init() {
    FixedBackOffPolicy smallerBackoff = new FixedBackOffPolicy();
    smallerBackoff.setBackOffPeriod(5L); // 5 ms
    template.setBackOffPolicy(smallerBackoff);
  }

  @Test
  void socketExceptionRetriesEventuallySucceed() throws Exception {
    MethodListResponse expectedResponse =
        new MethodListResponse().addMethodsItem(new MethodDetails().name("methodName"));

    CbasClient cbasClient = mock(CbasClient.class);
    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    doReturn(methodsApi).when(cbasClient).methodsApi(any(), any());

    assertEquals(expectedResponse, cbasService.getAllMethods(cbaseBaseUri, accessToken));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    MethodListResponse expectedResponse =
        new MethodListResponse().addMethodsItem(new MethodDetails().name("methodName"));

    CbasClient cbasClient = mock(CbasClient.class);
    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    doReturn(methodsApi).when(cbasClient).methodsApi(any(), any());

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          cbasService.getAllMethods(cbaseBaseUri, accessToken);
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    MethodListResponse expectedResponse =
        new MethodListResponse().addMethodsItem(new MethodDetails().name("methodName"));

    ApiException expectedException = new ApiException(400, "Bad Cbas");

    CbasClient cbasClient = mock(CbasClient.class);
    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(any(), any(), any()))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    doReturn(methodsApi).when(cbasClient).methodsApi(any(), any());

    CbasServiceApiException thrown =
        assertThrows(
            CbasServiceApiException.class,
            () -> {
              cbasService.getAllMethods(cbaseBaseUri, accessToken);
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void getAllMethodsTest() throws ApiException {
    MethodListResponse expectedResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name("methodName")
                    .methodId(UUID.randomUUID())
                    .description("this is my first method"))
            .addMethodsItem(
                new MethodDetails()
                    .name("secondMethodName")
                    .methodId(UUID.randomUUID())
                    .description("this is my second method"));

    CbasClient cbasClient = mock(CbasClient.class);
    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(any(), any(), any())).thenReturn(expectedResponse);

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    doReturn(methodsApi).when(cbasClient).methodsApi(any(), any());

    assertEquals(expectedResponse, cbasService.getAllMethods(cbaseBaseUri, accessToken));
  }

  @Test
  void createMethodTest() throws ApiException {
    PostMethodResponse expectedResponse = new PostMethodResponse().methodId(UUID.randomUUID());

    CbasClient cbasClient = mock(CbasClient.class);
    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.postMethod(any())).thenReturn(expectedResponse);

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    doReturn(methodsApi).when(cbasClient).methodsApi(any(), any());

    assertEquals(
        expectedResponse,
        cbasService.createMethod(cbaseBaseUri, accessToken, new PostMethodRequest()));
  }

  @Test
  void createRunSetTest() throws ApiException {
    RunSetStateResponse expectedResponse =
        new RunSetStateResponse()
            .state(RunSetState.RUNNING)
            .runSetId(UUID.randomUUID())
            .addRunsItem(new RunStateResponse().runId(UUID.randomUUID()).state(RunState.RUNNING));

    CbasClient cbasClient = mock(CbasClient.class);
    RunSetsApi runSetsApi = mock(RunSetsApi.class);
    when(runSetsApi.postRunSet(any())).thenReturn(expectedResponse);

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    doReturn(runSetsApi).when(cbasClient).runSetsApi(any(), any());

    assertEquals(
        expectedResponse, cbasService.createRunset(cbaseBaseUri, accessToken, new RunSetRequest()));
  }

  @Test
  void checkHealth() throws ApiException {
    CbasClient cbasClient = mock(CbasClient.class);
    PublicApi publicApi = mock(PublicApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    doReturn(publicApi).when(cbasClient).publicApi(any(), any());
    when(publicApi.getStatus()).thenReturn(systemStatus);

    CbasService cbasService = spy(new CbasService(cbasClient, template));
    HealthCheckWorkspaceApps.Result actualResult = cbasService.checkHealth("baseuri", "token");

    assertEquals(
        new HealthCheckWorkspaceApps.Result(systemStatus.isOk(), systemStatus.toString()),
        actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    CbasClient cbasClient = mock(CbasClient.class);
    PublicApi publicApi = mock(PublicApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);

    doReturn(publicApi).when(cbasClient).publicApi(any(), any());
    when(publicApi.getStatus()).thenThrow(apiException);

    HealthCheckWorkspaceApps.Result expectedResultOnFail =
        new HealthCheckWorkspaceApps.Result(false, apiException.getMessage());

    CbasService cbasService = spy(new CbasService(cbasClient, template));

    HealthCheckWorkspaceApps.Result actualResult = cbasService.checkHealth("baseuri", "token");

    assertEquals(expectedResultOnFail, actualResult);
  }
}
