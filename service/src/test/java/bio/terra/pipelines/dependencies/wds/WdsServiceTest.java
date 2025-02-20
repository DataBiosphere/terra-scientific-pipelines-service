package bio.terra.pipelines.dependencies.wds;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.DependencyNotAvailableException;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.net.SocketTimeoutException;
import java.util.List;
import org.databiosphere.workspacedata.api.GeneralWdsInformationApi;
import org.databiosphere.workspacedata.api.RecordsApi;
import org.databiosphere.workspacedata.api.SchemaApi;
import org.databiosphere.workspacedata.client.ApiException;
import org.databiosphere.workspacedata.model.RecordRequest;
import org.databiosphere.workspacedata.model.RecordResponse;
import org.databiosphere.workspacedata.model.RecordTypeSchema;
import org.databiosphere.workspacedata.model.StatusResponse;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.InjectMocks;
import org.mockito.junit.jupiter.MockitoExtension;
import org.mockito.stubbing.Answer;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

@ExtendWith(MockitoExtension.class)
class WdsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks WdsService wdsService;
  @MockitoBean WdsClient wdsClient;

  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

  private final String workspaceId = "workspaceId";
  private final String wdsApiVersion = "testapiv"; // matches value in test application.yaml
  private final String recordType = "recordType";
  private final String recordId = "recordId";
  private final String baseUri = "wdsUri";
  private final String bearerToken = "bearerToken";

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
    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(workspaceId, wdsApiVersion, recordType, recordId))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(recordsApi).when(wdsClient).recordsApi(baseUri, bearerToken);

    assertEquals(
        expectedResponse,
        wdsService.getRecord(baseUri, bearerToken, recordType, workspaceId, recordId));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(workspaceId, wdsApiVersion, recordType, recordId))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(recordsApi).when(wdsClient).recordsApi(baseUri, bearerToken);

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          wdsService.getRecord(baseUri, bearerToken, recordType, workspaceId, recordId);
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    ApiException expectedException = new ApiException(400, "Bad Wds");

    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(workspaceId, wdsApiVersion, recordType, recordId))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    doReturn(recordsApi).when(wdsClient).recordsApi(baseUri, bearerToken);

    WdsServiceApiException thrown =
        assertThrows(
            WdsServiceApiException.class,
            () -> {
              wdsService.getRecord(baseUri, bearerToken, recordType, workspaceId, recordId);
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void apiExceptionsDoNotRetryDependencyException() throws Exception {
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("typeTest");

    DependencyNotAvailableException expectedException =
        DependencyNotAvailableException.formatDependencyNotAvailableExceptionHelper(
            "WDS", "Bad Wds");

    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.getRecord(workspaceId, wdsApiVersion, recordType, recordId))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    doReturn(recordsApi).when(wdsClient).recordsApi(baseUri, bearerToken);

    WdsServiceNotAvailableException thrown =
        assertThrows(
            WdsServiceNotAvailableException.class,
            () -> {
              wdsService.getRecord(baseUri, bearerToken, recordType, workspaceId, recordId);
            });
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void updateRecord() throws Exception {
    RecordRequest recordRequest = new RecordRequest().putAdditionalProperty("foo", "bar");
    RecordResponse expectedResponse = new RecordResponse().id("idTest").type("updateRecordTest");

    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.updateRecord(recordRequest, workspaceId, wdsApiVersion, recordType, recordId))
        .thenReturn(expectedResponse);

    doReturn(recordsApi).when(wdsClient).recordsApi(baseUri, bearerToken);

    assertEquals(
        expectedResponse,
        wdsService.updateRecord(
            baseUri, bearerToken, recordRequest, workspaceId, recordType, recordId));
  }

  @Test
  void createOrReplaceRecord() throws Exception {
    RecordRequest recordRequest = new RecordRequest().putAdditionalProperty("foo", "baz");
    RecordResponse expectedResponse =
        new RecordResponse().id("idTest").type("createOrReplaceRecordTest");

    RecordsApi recordsApi = mock(RecordsApi.class);
    when(recordsApi.createOrReplaceRecord(
            recordRequest, workspaceId, wdsApiVersion, recordType, recordId, "primaryKey"))
        .thenReturn(expectedResponse);

    doReturn(recordsApi).when(wdsClient).recordsApi(baseUri, bearerToken);

    assertEquals(
        expectedResponse,
        wdsService.createOrReplaceRecord(
            baseUri, bearerToken, recordRequest, workspaceId, recordType, recordId, "primaryKey"));
  }

  @Test
  void querySchema() throws Exception {
    List<RecordTypeSchema> expectedResponse =
        List.of(
            new RecordTypeSchema().name("schemaTest"), new RecordTypeSchema().name("schemaTest2"));

    SchemaApi schemaApi = mock(SchemaApi.class);
    when(schemaApi.describeAllRecordTypes(workspaceId, wdsApiVersion)).thenReturn(expectedResponse);

    doReturn(schemaApi).when(wdsClient).schemaApi(baseUri, bearerToken);

    assertEquals(2, wdsService.querySchema(baseUri, bearerToken, workspaceId).size());
  }

  @Test
  void checkHealth() throws ApiException {
    GeneralWdsInformationApi generalWdsInformationApi = mock(GeneralWdsInformationApi.class);

    StatusResponse statusResponse = new StatusResponse().status("UP");

    when(wdsClient.generalWdsInformationApi(baseUri, bearerToken))
        .thenReturn(generalWdsInformationApi);
    when(generalWdsInformationApi.statusGet()).thenReturn(statusResponse);

    HealthCheckWorkspaceApps.Result actualResult = wdsService.checkHealth(baseUri, bearerToken);

    assertEquals(
        new HealthCheckWorkspaceApps.Result(true, statusResponse.toString()), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    GeneralWdsInformationApi generalWdsInformationApi = mock(GeneralWdsInformationApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(wdsClient.generalWdsInformationApi(baseUri, bearerToken))
        .thenReturn(generalWdsInformationApi);
    when(generalWdsInformationApi.statusGet()).thenThrow(apiException);

    HealthCheckWorkspaceApps.Result expectedResultOnFail =
        new HealthCheckWorkspaceApps.Result(false, apiException.getMessage());

    HealthCheckWorkspaceApps.Result actualResult = wdsService.checkHealth(baseUri, bearerToken);

    assertEquals(expectedResultOnFail, actualResult);
  }
}
