package bio.terra.pipelines.dependencies.rawls;

import static bio.terra.pipelines.testutils.TestUtils.VALID_METHOD_CONFIGURATION;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;
import static org.mockito.Mockito.doThrow;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.api.*;
import bio.terra.rawls.client.ApiException;
import bio.terra.rawls.model.*;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.net.SocketTimeoutException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
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
  @Autowired ObjectMapper objectMapper;

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

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

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

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

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

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

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
  void getWorkspaceDetails() throws Exception {
    String testBucketName = "bucketName";
    String testGoogleProject = "googleProject";
    WorkspaceDetails expectedWorkspaceDetails =
        new WorkspaceDetails().bucketName(testBucketName).googleProject(testGoogleProject);
    WorkspaceResponse expectedResponse =
        new WorkspaceResponse().workspace(expectedWorkspaceDetails);

    rawlsClient = mock(RawlsClient.class);
    WorkspacesApi workspacesApi = mock(WorkspacesApi.class);
    when(workspacesApi.listWorkspaceDetails(any(), any(), any())).thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

    doReturn(workspacesApi).when(rawlsClient).getWorkspacesApi(any());

    assertEquals(
        expectedWorkspaceDetails,
        rawlsService.getWorkspaceDetails("token", "workspaceNamespace", "workspace"));
  }

  @Test
  void getWorkspaceBucketName() {
    WorkspaceDetails workspaceDetails = new WorkspaceDetails().bucketName("bucketName");
    assertEquals("bucketName", rawlsService.getWorkspaceBucketName(workspaceDetails));
  }

  @Test
  void getWorkspaceBucketNameNullThrows() {
    WorkspaceDetails workspaceDetails = new WorkspaceDetails();
    assertThrows(
        InternalServerErrorException.class,
        () -> rawlsService.getWorkspaceBucketName(workspaceDetails));
  }

  @Test
  void getWorkspaceGoogleProject() {
    WorkspaceDetails workspaceDetails = new WorkspaceDetails().googleProject("googleProject");
    assertEquals("googleProject", rawlsService.getWorkspaceGoogleProject(workspaceDetails));
  }

  @Test
  void getWorkspaceGoogleProjectNullThrows() {
    WorkspaceDetails workspaceDetails = new WorkspaceDetails();
    assertThrows(
        InternalServerErrorException.class,
        () -> rawlsService.getWorkspaceGoogleProject(workspaceDetails));
  }

  @Test
  void submitWorkflow() throws Exception {
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    rawlsClient = mock(RawlsClient.class);
    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(any(), any(), any())).thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

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

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

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

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

    doReturn(entitiesApi).when(rawlsClient).getEntitiesApi(any());

    assertEquals(expectedResponse, rawlsService.upsertDataTableEntity(any(), any(), any(), any()));
  }

  @Test
  void getDataTableEntity() throws Exception {
    Entity expectedResponse = new Entity().name("entityName").entityType("entityType");

    rawlsClient = mock(RawlsClient.class);
    EntitiesApi entitiesApi = mock(EntitiesApi.class);
    when(entitiesApi.getEntity(any(), any(), any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

    doReturn(entitiesApi).when(rawlsClient).getEntitiesApi(any());

    assertEquals(
        expectedResponse,
        rawlsService.getDataTableEntity(
            "token", "billingProject", "workspace", "entityType", "entityName"));
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
  void getCurrentMethodConfig() throws Exception {
    MethodConfiguration expectedResponse = VALID_METHOD_CONFIGURATION;

    rawlsClient = mock(RawlsClient.class);
    MethodconfigsApi methodconfigsApi = mock(MethodconfigsApi.class);
    when(methodconfigsApi.getMethodConfiguration(any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

    doReturn(methodconfigsApi).when(rawlsClient).getMethodConfigsApi(any());

    assertEquals(
        expectedResponse, rawlsService.getCurrentMethodConfigForMethod(any(), any(), any(), any()));
  }

  @Test
  void updateMethodConfig() throws Exception {
    MethodConfiguration methodConfiguration = VALID_METHOD_CONFIGURATION;

    ValidatedMethodConfiguration expectedResponse =
        new ValidatedMethodConfiguration()
            .methodConfiguration(methodConfiguration)
            .invalidInputs(Collections.emptyMap());

    rawlsClient = mock(RawlsClient.class);
    MethodconfigsApi methodconfigsApi = mock(MethodconfigsApi.class);
    when(methodconfigsApi.updateMethodConfiguration(any(), any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    rawlsService = spy(new RawlsService(rawlsClient, template, objectMapper));

    doReturn(methodconfigsApi).when(rawlsClient).getMethodConfigsApi(any());

    assertEquals(
        expectedResponse, rawlsService.setMethodConfigForMethod(any(), any(), any(), any(), any()));
  }

  @Test
  void validateMethodConfig() {
    // generate values that should pass with VALID_METHOD_CONFIGURATION
    String workflowName = "workflowName";
    List<PipelineInputDefinition> expectedInputs =
        List.of(generatePipelineInputDefinitionWithWdlVariableName("first_input"));
    List<PipelineOutputDefinition> expectedOutputs =
        List.of(generatePipelineOutputDefinitionWithWdlVariableName("first_output"));
    String expectedMethodVersion = "1.2.3";

    assertTrue(
        rawlsService.validateMethodConfig(
            VALID_METHOD_CONFIGURATION,
            workflowName,
            expectedInputs,
            expectedOutputs,
            expectedMethodVersion));

    // different method version
    String differentMethodVersion = "1.1.1";
    assertFalse(
        rawlsService.validateMethodConfig(
            VALID_METHOD_CONFIGURATION,
            workflowName,
            expectedInputs,
            expectedOutputs,
            differentMethodVersion));

    // different expected inputs
    MethodConfiguration invalidInputMethodConfig =
        new MethodConfiguration()
            .name("name")
            .inputs(Map.of("workflowName.first_input", "this.wrong_input_reference"))
            .outputs(Map.of("workflowName.first_output", "this.first_output"))
            .methodRepoMethod(
                new MethodRepoMethod()
                    .methodName("methodName")
                    .methodNamespace("namespace")
                    .methodVersion("1.2.3")
                    .methodUri("this/is/a/uri/with/a/version/1.2.3"));
    // test wrong reference
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidInputMethodConfig,
            workflowName,
            List.of(generatePipelineInputDefinitionWithWdlVariableName("first_input")),
            expectedOutputs,
            expectedMethodVersion));
    // test missing key
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidInputMethodConfig,
            workflowName,
            List.of(generatePipelineInputDefinitionWithWdlVariableName("second_input")),
            expectedOutputs,
            expectedMethodVersion));

    // different expected outputs
    MethodConfiguration invalidOutputsMethodConfig =
        new MethodConfiguration()
            .name("name")
            .inputs(Map.of("workflowName.first_input", "this.first_input"))
            .outputs(Map.of("workflowName.first_output", "this.wrong_output_reference"))
            .methodRepoMethod(
                new MethodRepoMethod()
                    .methodName("methodName")
                    .methodNamespace("namespace")
                    .methodVersion("1.2.3")
                    .methodUri("this/is/a/uri/with/a/version/1.2.3"));
    // test wrong reference
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidOutputsMethodConfig,
            workflowName,
            expectedInputs,
            List.of(generatePipelineOutputDefinitionWithWdlVariableName("first_output")),
            expectedMethodVersion));
    // test missing key
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidOutputsMethodConfig,
            workflowName,
            expectedInputs,
            List.of(generatePipelineOutputDefinitionWithWdlVariableName("second_output")),
            expectedMethodVersion));
  }

  // take an "invalid" configuration and update it to match the expected method config
  @Test
  void updateMethodConfigToBeValid() {
    // generate values that should match the VALID_METHOD_CONFIGURATION
    String workflowName = "workflowName";
    List<PipelineInputDefinition> expectedInputs =
        List.of(generatePipelineInputDefinitionWithWdlVariableName("first_input"));
    List<PipelineOutputDefinition> expectedOutputs =
        List.of(generatePipelineOutputDefinitionWithWdlVariableName("first_output"));
    String expectedMethodVersion = "1.2.3";

    // generate method config that does not match the VALID_METHOD_CONFIGURATION
    MethodConfiguration invalidMethodConfig =
        new MethodConfiguration()
            .name("name")
            .inputs(Map.of("workflowName.first_input", "this.wrong_input_reference"))
            .outputs(Map.of("workflowName.first_output", "this.wrong_output_reference"))
            .methodRepoMethod(
                new MethodRepoMethod()
                    .methodName("methodName")
                    .methodNamespace("namespace")
                    .methodVersion("0.0.1")
                    .methodUri("this/is/a/uri/with/a/version/0.0.1"));

    // assert the two method configs are not equal initially
    assertNotEquals(VALID_METHOD_CONFIGURATION, invalidMethodConfig);

    MethodConfiguration updatedMethodConfig =
        rawlsService.updateMethodConfigToBeValid(
            invalidMethodConfig,
            workflowName,
            expectedInputs,
            expectedOutputs,
            expectedMethodVersion);
    assertEquals(VALID_METHOD_CONFIGURATION, updatedMethodConfig);
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

  private PipelineInputDefinition generatePipelineInputDefinitionWithWdlVariableName(
      String wdlVariableName) {
    return new PipelineInputDefinition(null, null, wdlVariableName, null, null, null, null, null);
  }

  private PipelineOutputDefinition generatePipelineOutputDefinitionWithWdlVariableName(
      String wdlVariableName) {
    return new PipelineOutputDefinition(null, null, wdlVariableName, null);
  }
}
