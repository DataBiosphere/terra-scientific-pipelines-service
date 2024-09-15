package bio.terra.pipelines.dependencies.workspacemanager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
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
import bio.terra.workspace.model.*;
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
    when(resourceApi.enumerateResources(
            workspaceId,
            null,
            null,
            ResourceType.AZURE_STORAGE_CONTAINER,
            StewardshipType.CONTROLLED))
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
    when(resourceApi.enumerateResources(
            workspaceId,
            null,
            null,
            ResourceType.AZURE_STORAGE_CONTAINER,
            StewardshipType.CONTROLLED))
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
    when(resourceApi.enumerateResources(
            workspaceId,
            null,
            null,
            ResourceType.AZURE_STORAGE_CONTAINER,
            StewardshipType.CONTROLLED))
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
    when(resourceApi.enumerateResources(
            workspaceId,
            null,
            null,
            ResourceType.AZURE_STORAGE_CONTAINER,
            StewardshipType.CONTROLLED))
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
    when(resourceApi.enumerateResources(
            workspaceId,
            null,
            null,
            ResourceType.AZURE_STORAGE_CONTAINER,
            StewardshipType.CONTROLLED))
        .thenReturn(expectedResponse);

    doReturn(resourceApi).when(workspaceManagerClient).getResourceApi(authToken);

    // we grab the resourceId from the first resource in the list
    assertEquals(
        resourceId, workspaceManagerService.getWorkspaceStorageResourceId(workspaceId, authToken));
  }

  @Test
  void getReadSasUrlForFile() throws Exception {
    UUID resourceId = UUID.randomUUID();

    String expectedSasUrl = "https://sas.url";
    CreatedAzureStorageContainerSasToken expectedSasToken =
        new CreatedAzureStorageContainerSasToken().url(expectedSasUrl);

    ControlledAzureResourceApi controlledAzureResourceApi = mock(ControlledAzureResourceApi.class);
    when(controlledAzureResourceApi.createAzureStorageContainerSasToken(
            workspaceId, resourceId, null, 86400L, "r", "some/blob/path"))
        .thenReturn(expectedSasToken);

    doReturn(controlledAzureResourceApi)
        .when(workspaceManagerClient)
        .getControlledAzureResourceApi(authToken);

    assertEquals(
        expectedSasUrl,
        workspaceManagerService.getReadSasUrlForBlob(
            workspaceId, resourceId, "some/blob/path", authToken));
  }

  @Test
  void getWriteSasUrlForFile() throws Exception {
    UUID resourceId = UUID.randomUUID();

    String expectedSasUrl = "https://sas.url";
    CreatedAzureStorageContainerSasToken expectedSasToken =
        new CreatedAzureStorageContainerSasToken().url(expectedSasUrl);

    ControlledAzureResourceApi controlledAzureResourceApi = mock(ControlledAzureResourceApi.class);
    when(controlledAzureResourceApi.createAzureStorageContainerSasToken(
            workspaceId, resourceId, null, 86400L, "w", "some/blob/path"))
        .thenReturn(expectedSasToken);

    doReturn(controlledAzureResourceApi)
        .when(workspaceManagerClient)
        .getControlledAzureResourceApi(authToken);

    assertEquals(
        expectedSasUrl,
        workspaceManagerService.getWriteSasUrlForBlob(
            workspaceId, resourceId, "some/blob/path", authToken));
  }
}
