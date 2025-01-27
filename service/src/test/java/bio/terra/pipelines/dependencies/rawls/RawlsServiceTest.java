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

    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(null, "workspaceNamespace", "workspaceName"))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi("accessToken");

    assertEquals(
        expectedResponse,
        rawlsService.submitWorkflow("accessToken", null, "workspaceNamespace", "workspaceName"));
  }

  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(null, "workspaceNamespace", "workspaceName"))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi("accessToken");

    assertThrows(
        SocketTimeoutException.class,
        () -> {
          rawlsService.submitWorkflow("accessToken", null, "workspaceNamespace", "workspaceName");
        });
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    SubmissionReport expectedResponse =
        new SubmissionReport().status("status").submissionId(UUID.randomUUID().toString());

    ApiException expectedException = new ApiException(400, "Bad Rawls");
    SubmissionRequest emptySubmissionRequest = new SubmissionRequest();

    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(emptySubmissionRequest, "workspaceNamespace", "workspace"))
        .thenThrow(expectedException)
        .thenReturn(expectedResponse);

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi("blah");
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

    WorkspacesApi workspacesApi = mock(WorkspacesApi.class);
    when(workspacesApi.listWorkspaceDetails(eq("workspaceNamespace"), eq("workspace"), anyList()))
        .thenReturn(expectedResponse);

    doReturn(workspacesApi).when(rawlsClient).getWorkspacesApi("token");

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

    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.createSubmission(null, "workspaceNamespace", "workspaceName"))
        .thenReturn(expectedResponse);

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi("accessToken");

    assertEquals(
        expectedResponse,
        rawlsService.submitWorkflow("accessToken", null, "workspaceNamespace", "workspaceName"));
  }

  @Test
  void getSubmissionStatus() throws Exception {
    UUID expectedUUID = UUID.randomUUID();
    Submission expectedResponse =
        new Submission()
            .status(SubmissionStatus.SUBMITTING)
            .submissionId(expectedUUID.toString())
            .cost(3.45F);

    SubmissionsApi submissionsApi = mock(SubmissionsApi.class);
    when(submissionsApi.getSubmissionStatus(
            "workspaceNamespace", "workspace", expectedUUID.toString()))
        .thenReturn(expectedResponse);

    doReturn(submissionsApi).when(rawlsClient).getSubmissionsApi("token");

    assertEquals(
        expectedResponse,
        rawlsService.getSubmissionStatus("token", "workspaceNamespace", "workspace", expectedUUID));
  }

  @Test
  void upsertDataTableEntity() throws Exception {
    Entity expectedResponse = new Entity().name("entityName").entityType("entityType");

    EntitiesApi entitiesApi = mock(EntitiesApi.class);
    when(entitiesApi.createEntity(null, "workspaceNamespace", "workspaceName"))
        .thenReturn(expectedResponse);

    doReturn(entitiesApi).when(rawlsClient).getEntitiesApi("accessToken");

    assertEquals(
        expectedResponse,
        rawlsService.upsertDataTableEntity(
            "accessToken", "workspaceNamespace", "workspaceName", null));
  }

  @Test
  void getDataTableEntity() throws Exception {
    Entity expectedResponse = new Entity().name("entityName").entityType("entityType");

    EntitiesApi entitiesApi = mock(EntitiesApi.class);
    when(entitiesApi.getEntity(
            "billingProject", "workspace", "entityType", "entityName", null, null))
        .thenReturn(expectedResponse);

    doReturn(entitiesApi).when(rawlsClient).getEntitiesApi("token");

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

    MethodconfigsApi methodconfigsApi = mock(MethodconfigsApi.class);
    when(methodconfigsApi.getMethodConfiguration(
            "workspaceNamespace", "workspaceName", "workspaceNamespace", "methodName"))
        .thenReturn(expectedResponse);

    doReturn(methodconfigsApi).when(rawlsClient).getMethodConfigsApi("accessToken");

    assertEquals(
        expectedResponse,
        rawlsService.getCurrentMethodConfigForMethod(
            "accessToken", "workspaceNamespace", "workspaceName", "methodName"));
  }

  @Test
  void updateMethodConfig() throws Exception {
    ValidatedMethodConfiguration expectedResponse =
        new ValidatedMethodConfiguration()
            .methodConfiguration(VALID_METHOD_CONFIGURATION)
            .invalidInputs(Collections.emptyMap());

    MethodconfigsApi methodconfigsApi = mock(MethodconfigsApi.class);
    when(methodconfigsApi.updateMethodConfiguration(
            null, "workspaceNamespace", "workspaceName", "workspaceNamespace", "methodName"))
        .thenReturn(expectedResponse);

    doReturn(methodconfigsApi).when(rawlsClient).getMethodConfigsApi("accessToken");

    assertEquals(
        expectedResponse,
        rawlsService.setMethodConfigForMethod(
            "accessToken", null, "workspaceNamespace", "workspaceName", "methodName"));
  }

  @Test
  void validateMethodConfig() {
    // generate values that should pass with VALID_METHOD_CONFIGURATION
    String expectedWorkflowName = "workflowName";
    String expectedDataTableEntity = "imputation_beagle";
    List<PipelineInputDefinition> expectedInputs =
        List.of(generatePipelineInputDefinitionWithWdlVariableName("first_input"));
    List<PipelineOutputDefinition> expectedOutputs =
        List.of(generatePipelineOutputDefinitionWithWdlVariableName("first_output"));
    String expectedMethodVersion = "1.2.3";

    assertTrue(
        rawlsService.validateMethodConfig(
            VALID_METHOD_CONFIGURATION,
            expectedDataTableEntity,
            expectedWorkflowName,
            expectedInputs,
            expectedOutputs,
            expectedMethodVersion));

    // different data table entity name
    String differentDataTableEntityName = "different_table";
    assertFalse(
        rawlsService.validateMethodConfig(
            VALID_METHOD_CONFIGURATION,
            differentDataTableEntityName,
            expectedWorkflowName,
            expectedInputs,
            expectedOutputs,
            expectedMethodVersion));

    // different method version
    String differentMethodVersion = "1.1.1";
    assertFalse(
        rawlsService.validateMethodConfig(
            VALID_METHOD_CONFIGURATION,
            expectedDataTableEntity,
            expectedWorkflowName,
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
                    .methodUri("this/is/a/uri/with/a/version/1.2.3"))
            .rootEntityType("imputation_beagle");
    // test wrong reference
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidInputMethodConfig,
            expectedDataTableEntity,
            expectedWorkflowName,
            List.of(generatePipelineInputDefinitionWithWdlVariableName("first_input")),
            expectedOutputs,
            expectedMethodVersion));
    // test missing key
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidInputMethodConfig,
            expectedDataTableEntity,
            expectedWorkflowName,
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
                    .methodUri("this/is/a/uri/with/a/version/1.2.3"))
            .rootEntityType("imputation_beagle");
    // test wrong reference
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidOutputsMethodConfig,
            expectedDataTableEntity,
            expectedWorkflowName,
            expectedInputs,
            List.of(generatePipelineOutputDefinitionWithWdlVariableName("first_output")),
            expectedMethodVersion));
    // test missing key
    assertFalse(
        rawlsService.validateMethodConfig(
            invalidOutputsMethodConfig,
            expectedDataTableEntity,
            expectedWorkflowName,
            expectedInputs,
            List.of(generatePipelineOutputDefinitionWithWdlVariableName("second_output")),
            expectedMethodVersion));
  }

  // take an "invalid" configuration and update it to match the expected method config
  @Test
  void updateMethodConfigToBeValid() {
    // generate values that should match the VALID_METHOD_CONFIGURATION
    String expectedWorkflowName = "workflowName";
    String expectedDataTableEntity = "imputation_beagle";
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
                    .methodUri("this/is/a/uri/with/a/version/0.0.1"))
            .rootEntityType("different_from_valid");

    // assert the two method configs are not equal initially
    assertNotEquals(VALID_METHOD_CONFIGURATION, invalidMethodConfig);

    MethodConfiguration updatedMethodConfig =
        rawlsService.updateMethodConfigToBeValid(
            invalidMethodConfig,
            expectedDataTableEntity,
            expectedWorkflowName,
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
    return new PipelineInputDefinition(
        null, null, wdlVariableName, null, null, true, true, false, null);
  }

  private PipelineOutputDefinition generatePipelineOutputDefinitionWithWdlVariableName(
      String wdlVariableName) {
    return new PipelineOutputDefinition(null, null, wdlVariableName, null);
  }
}
