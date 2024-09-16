package bio.terra.pipelines.dependencies.cbas;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;
import static org.mockito.Mockito.doReturn;

import bio.terra.cbas.api.MethodsApi;
import bio.terra.cbas.api.PublicApi;
import bio.terra.cbas.api.RunSetsApi;
import bio.terra.cbas.api.RunsApi;
import bio.terra.cbas.client.ApiException;
import bio.terra.cbas.model.*;
import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.net.SocketTimeoutException;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.InjectMocks;
import org.mockito.junit.jupiter.MockitoExtension;
import org.mockito.stubbing.Answer;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

@ExtendWith(MockitoExtension.class)
class CbasServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks CbasService cbasService;
  @MockBean CbasClient cbasClient;

  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  final CbasConfiguration cbasConfiguration = new CbasConfiguration();

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

    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(null, null, null))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(methodsApi).when(cbasClient).methodsApi(cbaseBaseUri, accessToken);

    assertEquals(expectedResponse, cbasService.getAllMethods(cbaseBaseUri, accessToken));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    MethodListResponse expectedResponse =
        new MethodListResponse().addMethodsItem(new MethodDetails().name("methodName"));

    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(null, null, null))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(methodsApi).when(cbasClient).methodsApi(cbaseBaseUri, accessToken);

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

    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(null, null, null))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    doReturn(methodsApi).when(cbasClient).methodsApi(cbaseBaseUri, accessToken);

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

    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.getMethods(null, null, null)).thenReturn(expectedResponse);

    doReturn(methodsApi).when(cbasClient).methodsApi(cbaseBaseUri, accessToken);

    assertEquals(expectedResponse, cbasService.getAllMethods(cbaseBaseUri, accessToken));
  }

  @Test
  void createMethodTest() throws ApiException {
    PostMethodResponse expectedResponse = new PostMethodResponse().methodId(UUID.randomUUID());
    PostMethodRequest postMethodRequest = new PostMethodRequest().methodName("hi name");
    MethodsApi methodsApi = mock(MethodsApi.class);
    when(methodsApi.postMethod(postMethodRequest)).thenReturn(expectedResponse);

    doReturn(methodsApi).when(cbasClient).methodsApi(cbaseBaseUri, accessToken);

    assertEquals(
        expectedResponse, cbasService.createMethod(cbaseBaseUri, accessToken, postMethodRequest));
  }

  @Test
  void createRunSetTest() throws ApiException {
    RunSetStateResponse expectedResponse =
        new RunSetStateResponse()
            .state(RunSetState.RUNNING)
            .runSetId(UUID.randomUUID())
            .addRunsItem(new RunStateResponse().runId(UUID.randomUUID()).state(RunState.RUNNING));
    RunSetRequest runSetRequest = new RunSetRequest().runSetName("hi run ste name");

    RunSetsApi runSetsApi = mock(RunSetsApi.class);
    when(runSetsApi.postRunSet(runSetRequest)).thenReturn(expectedResponse);

    doReturn(runSetsApi).when(cbasClient).runSetsApi(cbaseBaseUri, accessToken);

    assertEquals(
        expectedResponse, cbasService.createRunSet(cbaseBaseUri, accessToken, runSetRequest));
  }

  @Test
  void getRunsForRunSetTest() throws ApiException {
    RunLogResponse expectedResponse =
        new RunLogResponse().addRunsItem(new RunLog().state(RunState.RUNNING));
    UUID runSetId = UUID.randomUUID();

    RunsApi runsApi = mock(RunsApi.class);
    when(runsApi.getRuns(runSetId)).thenReturn(expectedResponse);

    doReturn(runsApi).when(cbasClient).runsApi(cbaseBaseUri, accessToken);

    assertEquals(
        expectedResponse, cbasService.getRunsForRunSet(cbaseBaseUri, accessToken, runSetId));
  }

  @Test
  void checkHealth() throws ApiException {
    PublicApi publicApi = mock(PublicApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    doReturn(publicApi).when(cbasClient).publicApi(cbaseBaseUri, accessToken);
    when(publicApi.getStatus()).thenReturn(systemStatus);

    HealthCheckWorkspaceApps.Result actualResult =
        cbasService.checkHealth(cbaseBaseUri, accessToken);

    assertEquals(
        new HealthCheckWorkspaceApps.Result(systemStatus.isOk(), systemStatus.toString()),
        actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    PublicApi publicApi = mock(PublicApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);

    doReturn(publicApi).when(cbasClient).publicApi(cbaseBaseUri, accessToken);
    when(publicApi.getStatus()).thenThrow(apiException);

    HealthCheckWorkspaceApps.Result expectedResultOnFail =
        new HealthCheckWorkspaceApps.Result(false, apiException.getMessage());

    HealthCheckWorkspaceApps.Result actualResult =
        cbasService.checkHealth(cbaseBaseUri, accessToken);

    assertEquals(expectedResultOnFail, actualResult);
  }

  @Test
  void getMethodVersionIdFromMethodListResponse() {
    UUID methodVersionIdActual = UUID.randomUUID();
    MethodListResponse getAllMethodsResponse =
        new MethodListResponse()
            .addMethodsItem(
                new MethodDetails()
                    .name("random name that doesnt match anything")
                    .addMethodVersionsItem(
                        new MethodVersionDetails().methodVersionId(methodVersionIdActual)));
    UUID methodVersionIdNoMatch =
        CbasService.getMethodVersionIdFromMethodListResponse(
            getAllMethodsResponse, "notRandomName");
    assertNull(methodVersionIdNoMatch);

    UUID methodVersionIdExpected =
        CbasService.getMethodVersionIdFromMethodListResponse(
            getAllMethodsResponse, "random name that doesnt match anything");
    assertEquals(methodVersionIdActual, methodVersionIdExpected);
  }

  @Test
  void containsRunningLog() {
    RunLogResponse runningLog =
        new RunLogResponse()
            .addRunsItem(new RunLog().state(RunState.RUNNING))
            .addRunsItem(new RunLog().state(RunState.COMPLETE));
    assertTrue(CbasService.containsRunningRunLog(runningLog));

    RunLogResponse noRunningLog =
        new RunLogResponse()
            .addRunsItem(new RunLog().state(RunState.EXECUTOR_ERROR))
            .addRunsItem(new RunLog().state(RunState.COMPLETE));
    assertFalse(CbasService.containsRunningRunLog(noRunningLog));
  }
}
