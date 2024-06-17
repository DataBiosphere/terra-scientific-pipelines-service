package bio.terra.pipelines.dependencies.workspacemanager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.doNothing;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.workspace.api.ControlledAzureResourceApi;
import bio.terra.workspace.api.ResourceApi;
import bio.terra.workspace.api.UnauthenticatedApi;
import bio.terra.workspace.client.ApiException;
import bio.terra.workspace.model.CreatedAzureStorageContainerSasToken;
import bio.terra.workspace.model.ResourceDescription;
import bio.terra.workspace.model.ResourceList;
import bio.terra.workspace.model.ResourceMetadata;
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

class WorkspaceManagerServiceTest extends BaseEmbeddedDbTest {

  @Autowired @InjectMocks WorkspaceManagerService workspaceManagerService;
  @MockBean WorkspaceManagerClient workspaceManagerClient;

  final UUID workspaceId = UUID.randomUUID();

  final RetryConfiguration retryConfig = new RetryConfiguration();
  RetryTemplate template = retryConfig.listenerResetRetryTemplate();

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
  void checkHealth() throws ApiException {
    UnauthenticatedApi unauthenticatedApi = mock(UnauthenticatedApi.class);

    when(workspaceManagerClient.getUnauthenticatedApi()).thenReturn(unauthenticatedApi);
    doNothing()
        .when(unauthenticatedApi)
        .serviceStatus(); // Workspace Manager's serviceStatus() is a void method

    HealthCheck.Result actualResult = workspaceManagerService.checkHealth();

    assertEquals(new HealthCheck.Result(true, "Workspace Manager is ok"), actualResult);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    UnauthenticatedApi unauthenticatedApi = mock(UnauthenticatedApi.class);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);

    doReturn(unauthenticatedApi).when(workspaceManagerClient).getUnauthenticatedApi();
    doThrow(apiException).when(unauthenticatedApi).serviceStatus();

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());

    HealthCheck.Result actualResult = workspaceManagerService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);
  }

  @Test
  void socketExceptionRetriesEventuallySucceed() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));

    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);

    assertEquals(
        resourceId, workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));

    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);

    assertThrows(
        SocketTimeoutException.class,
        () -> workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    ApiException expectedException = new ApiException(400, "Bad Workspace manager");

    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenThrow(expectedException);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);

    WorkspaceManagerServiceApiException thrown =
        assertThrows(
            WorkspaceManagerServiceApiException.class,
            () -> workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, authToken));
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void getWorkspaceStorageResourceId() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));

    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);

    assertEquals(
        resourceId, workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  @Test
  void getWorkspaceStorageResourceIdMultipleWorkspaceStorageResources() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)))
            .addResourcesItem(
                new ResourceDescription()
                    .metadata(new ResourceMetadata().resourceId(UUID.randomUUID())));

    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);

    // we grab the resourceId from the first resource in the list
    assertEquals(
        resourceId, workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  @Test
  void getSasTokenForFile() throws Exception {
    UUID resourceId = UUID.randomUUID();
    String filePathFromWorkspace =
        "https://lzsomething.blob.core.windows.net/sc-%s/some/file/path".formatted(workspaceId);
    ResourceList expectedResourceListResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));
    String expectedSasUrl = "https://sas.url";
    CreatedAzureStorageContainerSasToken expectedSasToken =
        new CreatedAzureStorageContainerSasToken().url(expectedSasUrl);

    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenReturn(expectedResourceListResponse);

    ControlledAzureResourceApi controlledAzureResourceApi = mock(ControlledAzureResourceApi.class);
    when(controlledAzureResourceApi.createAzureStorageContainerSasToken(
            eq(workspaceId), eq(resourceId), any(), any(), any(), any()))
        .thenReturn(expectedSasToken);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);
    doReturn(controlledAzureResourceApi)
        .when(workspaceManagerClient)
        .getControlledAzureResourceApi(authToken);

    assertEquals(
        expectedSasUrl,
        workspaceManagerService.getSasTokenForFile(
            workspaceId, filePathFromWorkspace, "w", authToken));
  }

  @Test
  void getBlobNameFromFullPath() {
    String fullPath =
        "https://lze96253b07f13c61ef712bb.blob.core.windows.net/sc-68a43bd8-e744-4f1e-87a5-c44ecef157a3/workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    UUID controlWorkspaceId = UUID.fromString("68a43bd8-e744-4f1e-87a5-c44ecef157a3");
    String expectedBlobName =
        "workspace-services/cbas/terra-app-b1740821-d6e9-44b5-b53b-960953dea218/ImputationBeagle/1adb690d-3d02-4d4a-9dfa-17a31edd74f3/call-WriteEmptyFile/cacheCopy/execution/empty_file";
    assertEquals(
        expectedBlobName,
        workspaceManagerService.getBlobNameFromHttpUrl(fullPath, controlWorkspaceId));
  }
}
