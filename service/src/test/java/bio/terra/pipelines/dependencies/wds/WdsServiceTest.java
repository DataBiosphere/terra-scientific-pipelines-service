package bio.terra.pipelines.dependencies.wds;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.app.configuration.external.WdsServerConfiguration;
import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.DependencyNotAvailableException;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import java.net.SocketTimeoutException;
import java.util.List;
import org.databiosphere.workspacedata.api.GeneralWdsInformationApi;
import org.databiosphere.workspacedata.api.RecordsApi;
import org.databiosphere.workspacedata.api.SchemaApi;
import org.databiosphere.workspacedata.client.ApiException;
import org.databiosphere.workspacedata.model.RecordResponse;
import org.databiosphere.workspacedata.model.RecordTypeSchema;
import org.databiosphere.workspacedata.model.StatusResponse;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.junit.jupiter.MockitoExtension;
import org.mockito.stubbing.Answer;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

@ExtendWith(MockitoExtension.class)
class WdsServiceTest {
  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  final WdsServerConfiguration wdsServerConfiguration = new WdsServerConfiguration("v0.2", false);

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
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    WdsClient wdsClient = mock(WdsClient.class);
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(any(), any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(recordsApi).when(wdsClient).recordsApi(any(), any());

    assertEquals(expectedResponse, wdsService.getRecord(any(), any(), any(), any(), any()));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    WdsClient wdsClient = mock(WdsClient.class);
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(any(), any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(recordsApi).when(wdsClient).recordsApi(any(), any());

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          wdsService.getRecord(any(), any(), any(), any(), any());
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    ApiException expectedException = new ApiException(400, "Bad Wds");

    WdsClient wdsClient = mock(WdsClient.class);
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(any(), any(), any(), any()))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(recordsApi).when(wdsClient).recordsApi(any(), any());

    WdsServiceApiException thrown =
        assertThrows(
            WdsServiceApiException.class,
            () -> {
              wdsService.getRecord(any(), any(), any(), any(), any());
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void apiExceptionsDoNotRetryDependencyException() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    DependencyNotAvailableException expectedException =
        DependencyNotAvailableException.formatDependencyNotAvailableExceptionHelper(
            "WDS", "Bad Wds");

    WdsClient wdsClient = mock(WdsClient.class);
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(any(), any(), any(), any()))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(recordsApi).when(wdsClient).recordsApi(any(), any());

    WdsServiceNotAvailableException thrown =
        assertThrows(
            WdsServiceNotAvailableException.class,
            () -> {
              wdsService.getRecord(any(), any(), any(), any(), any());
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void updateRecord() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("updateRecordTest");

    WdsClient wdsClient = mock(WdsClient.class);
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.updateRecord(any(), any(), any(), any(), any())).thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(recordsApi).when(wdsClient).recordsApi(any(), any());

    assertEquals(
        expectedResponse, wdsService.updateRecord(any(), any(), any(), any(), any(), any()));
  }

  @Test
  void createOrReplaceRecord() throws Exception {
    RecordResponse expectedResponse =
        new RecordResponse().id("idTest").type("createOrReplaceRecordTest");

    WdsClient wdsClient = mock(WdsClient.class);
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.createOrReplaceRecord(any(), any(), any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(recordsApi).when(wdsClient).recordsApi(any(), any());

    assertEquals(
        expectedResponse,
        wdsService.createOrReplaceRecord(any(), any(), any(), any(), any(), any(), any()));
  }

  @Test
  void querySchema() throws Exception {
    List<RecordTypeSchema> expectedResponse =
        List.of(
            new RecordTypeSchema().name("schemaTest"), new RecordTypeSchema().name("schemaTest2"));

    WdsClient wdsClient = mock(WdsClient.class);
    SchemaApi schemaApi = mock(SchemaApi.class);
    when(schemaApi.describeAllRecordTypes(any(), any())).thenReturn(expectedResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    doReturn(schemaApi).when(wdsClient).schemaApi(any(), any());

    assertEquals(2, wdsService.querySchema(any(), any(), any()).size());
  }

  @Test
  void checkHealth() throws ApiException {
    WdsClient wdsClient = mock(WdsClient.class);
    GeneralWdsInformationApi generalWdsInformationApi = mock(GeneralWdsInformationApi.class);

    StatusResponse statusResponse = new StatusResponse().status("UP");

    when(wdsClient.generalWdsInformationApi(any(), any())).thenReturn(generalWdsInformationApi);
    when(generalWdsInformationApi.statusGet()).thenReturn(statusResponse);

    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    HealthCheckWorkspaceApps.Result actualResult = wdsService.checkHealth("baseuri", "token");

    assertEquals(
        new HealthCheckWorkspaceApps.Result(true, statusResponse.toString()), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    WdsClient wdsClient = mock(WdsClient.class);
    GeneralWdsInformationApi generalWdsInformationApi = mock(GeneralWdsInformationApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(wdsClient.generalWdsInformationApi(any(), any())).thenReturn(generalWdsInformationApi);
    when(generalWdsInformationApi.statusGet()).thenThrow(apiException);

    HealthCheckWorkspaceApps.Result expectedResultOnFail =
        new HealthCheckWorkspaceApps.Result(false, apiException.getMessage());
    WdsService wdsService = spy(new WdsService(wdsClient, wdsServerConfiguration, template));

    HealthCheckWorkspaceApps.Result actualResult = wdsService.checkHealth("baseuri", "token");

    assertEquals(expectedResultOnFail, actualResult);
  }
}
