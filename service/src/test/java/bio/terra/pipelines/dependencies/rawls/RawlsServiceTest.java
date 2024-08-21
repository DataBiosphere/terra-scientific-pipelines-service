package bio.terra.pipelines.dependencies.rawls;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;
import static org.mockito.Mockito.doThrow;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.api.EntitiesApi;
import bio.terra.rawls.api.StatusApi;
import bio.terra.rawls.api.SubmissionsApi;
import bio.terra.rawls.client.ApiException;
import bio.terra.rawls.model.*;
import java.net.SocketTimeoutException;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.stubbing.Answer;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

class RawlsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks RawlsService rawlsService;
  @MockBean RawlsClient rawlsClient;

  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

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
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    rawlsClient = mock(RawlsClient.class);
    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template));

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi(any());

    assertEquals(expectedResponse, rawlsService.submitWorkflow(any(), any(), any(), any()));
  }

  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    rawlsClient = mock(RawlsClient.class);
    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template));

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi(any());

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          rawlsService.submitWorkflow(any(), any(), any(), any());
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    ApiException expectedException = new ApiException(400, "Bad Rawls");

    rawlsClient = mock(RawlsClient.class);
    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(any(), any(), any()))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template));

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi(any());
    SubmissionRequest emptySubmissionRequest = new SubmissionRequest();
    RawlsServiceApiException thrown =
        assertThrows(
            RawlsServiceApiException.class,
            () ->
                rawlsService.submitWorkflow(
                    "blah", emptySubmissionRequest, "workspaceNamespace", "workspace"));
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void submitWorkflow() throws Exception {
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    rawlsClient = mock(RawlsClient.class);
    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(any(), any(), any())).thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template));

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi(any());

    assertEquals(expectedResponse, rawlsService.submitWorkflow(any(), any(), any(), any()));
  }

  @Test
  void getSubmissionStatus() throws Exception {
    UUID expectedUUID = UUID.randomUUID();
    Submission expectedResponse =
        new Submission()
            .status(SubmissionStatus.SUBMITTING)
            .submissionId(expectedUUID.toString())
            .cost(3.45F);

    rawlsClient = mock(RawlsClient.class);
    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.getSubmissionStatus(any(), any(), eq(expectedUUID.toString())))
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template));

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi(any());

    assertEquals(
        expectedResponse,
        rawlsService.getSubmissionStatus("token", "workspaceNamespace", "workspace", expectedUUID));
  }

  @Test
  void upsertDataTableEntity() throws Exception {
    Entity expectedResponse = new Entity().name("entityName").entityType("entityType");

    rawlsClient = mock(RawlsClient.class);
    EntitiesApi entitiesApi = mock(EntitiesApi.class);
    when(entitiesApi.createEntity(any(), any(), any())).thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template));

    doReturn(entitiesApi).when(rawlsClient).getEntitiesApi(any());

    assertEquals(expectedResponse, rawlsService.upsertDataTableEntity(any(), any(), any(), any()));
  }

  @Test
  void submissionIsNotRunning() {
    assertFalse(RawlsService.submissionIsRunning(new Submission().status(SubmissionStatus.DONE)));
    assertFalse(
        RawlsService.submissionIsRunning(new Submission().status(SubmissionStatus.ABORTED)));
  }

  @Test
  void submissionIsRunning() {
    assertTrue(
        RawlsService.submissionIsRunning(new Submission().status(SubmissionStatus.SUBMITTING)));
    assertTrue(
        RawlsService.submissionIsRunning(new Submission().status(SubmissionStatus.ABORTING)));
    assertTrue(
        RawlsService.submissionIsRunning(new Submission().status(SubmissionStatus.EVALUATING)));
    assertTrue(
        RawlsService.submissionIsRunning(new Submission().status(SubmissionStatus.SUBMITTED)));
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
