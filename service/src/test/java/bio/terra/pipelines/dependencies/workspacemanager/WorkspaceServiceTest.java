package bio.terra.pipelines.dependencies.workspacemanager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.spy;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import bio.terra.workspace.api.ControlledAzureResourceApi;
import bio.terra.workspace.api.ResourceApi;
import bio.terra.workspace.client.ApiException;
import bio.terra.workspace.model.CreatedAzureStorageContainerSasToken;
import bio.terra.workspace.model.ResourceDescription;
import bio.terra.workspace.model.ResourceList;
import bio.terra.workspace.model.ResourceMetadata;
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
public class WorkspaceServiceTest {

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
  void socketExceptionRetriesEventuallySucceed() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));

    WorkspaceClient workspaceClient = mock(WorkspaceClient.class);
    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    WorkspaceService workspaceService = spy(new WorkspaceService(workspaceClient, template));
    doReturn(resourceApi).when(workspaceService).getResourceApi(authToken);

    assertEquals(
        resourceId, workspaceService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  // our retry template only attempts a retryable call 3 total times
  @Test
  void socketExceptionRetriesEventuallyFail() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));

    WorkspaceClient workspaceClient = mock(WorkspaceClient.class);
    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenAnswer(errorAnswer)
        .thenReturn(expectedResponse);

    WorkspaceService workspaceService = spy(new WorkspaceService(workspaceClient, template));
    doReturn(resourceApi).when(workspaceService).getResourceApi(authToken);

    assertThrows(
        SocketTimeoutException.class,
        () -> workspaceService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  @Test
  void apiExceptionsDoNotRetry() throws Exception {
    ApiException expectedException = new ApiException(400, "Bad Workspace manager");

    WorkspaceClient workspaceClient = mock(WorkspaceClient.class);
    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenThrow(expectedException);

    WorkspaceService workspaceService = spy(new WorkspaceService(workspaceClient, template));
    doReturn(resourceApi).when(workspaceService).getResourceApi(authToken);

    WorkspaceServiceApiException thrown =
        assertThrows(
            WorkspaceServiceApiException.class,
            () -> workspaceService.getWorkspaceStorageResourceId(workspaceId, authToken));
    assertEquals(expectedException, thrown.getCause());
  }

  @Test
  void getWorkspaceStorageResourceId() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));

    WorkspaceClient workspaceClient = mock(WorkspaceClient.class);
    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenReturn(expectedResponse);

    WorkspaceService workspaceService = spy(new WorkspaceService(workspaceClient, template));
    doReturn(resourceApi).when(workspaceService).getResourceApi(authToken);

    assertEquals(
        resourceId, workspaceService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  @Test
  void getSasTokenForFile() throws Exception {
    UUID resourceId = UUID.randomUUID();
    ResourceList expectedResourceListResponse =
        new ResourceList()
            .addResourcesItem(
                new ResourceDescription().metadata(new ResourceMetadata().resourceId(resourceId)));
    String expectedSasUrl = "https://sas.url";
    CreatedAzureStorageContainerSasToken expectedSasToken =
        new CreatedAzureStorageContainerSasToken().url(expectedSasUrl);

    WorkspaceClient workspaceClient = mock(WorkspaceClient.class);
    ResourceApi resourceApi = mock(ResourceApi.class);
    when(resourceApi.enumerateResources(eq(workspaceId), any(), any(), any(), any()))
        .thenReturn(expectedResourceListResponse);

    ControlledAzureResourceApi controlledAzureResourceApi = mock(ControlledAzureResourceApi.class);
    when(controlledAzureResourceApi.createAzureStorageContainerSasToken(
            eq(workspaceId), eq(resourceId), any(), any(), any(), any()))
        .thenReturn(expectedSasToken);

    WorkspaceService workspaceService = spy(new WorkspaceService(workspaceClient, template));
    doReturn(resourceApi).when(workspaceService).getResourceApi(authToken);
    doReturn(controlledAzureResourceApi)
        .when(workspaceService)
        .getControlledAzureResourceApi(authToken);

    assertEquals(
        expectedSasUrl,
        workspaceService.getSasTokenForFile(workspaceId, "blobName", "w", authToken));
  }
}
