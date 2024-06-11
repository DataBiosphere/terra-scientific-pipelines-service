package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import bio.terra.workspace.api.ControlledAzureResourceApi;
import bio.terra.workspace.api.ResourceApi;
import bio.terra.workspace.api.UnauthenticatedApi;
import bio.terra.workspace.client.ApiException;
import bio.terra.workspace.model.CreatedAzureStorageContainerSasToken;
import bio.terra.workspace.model.ResourceList;
import bio.terra.workspace.model.ResourceType;
import bio.terra.workspace.model.StewardshipType;
import java.util.UUID;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class WorkspaceService implements HealthCheck {
  private final WorkspaceClient workspaceClient;
  private final RetryTemplate listenerResetRetryTemplate;

  public WorkspaceService(
      WorkspaceClient workspaceClient, RetryTemplate listenerResetRetryTemplate) {
    this.workspaceClient = workspaceClient;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  public Result checkHealth() {
    // No access token needed since this is an unauthenticated API.
    try {
      getUnauthenticatedApi()
          .serviceStatus(); // Workspace Manager's serviceStatus() is a void method
      return new Result(true, "Workspace Manager is ok");
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  public ApiSystemStatusSystems checkHealthApiSystemStatus() {
    Result healthResult = checkHealth();
    return new ApiSystemStatusSystems()
        .ok(healthResult.isOk())
        .addMessagesItem(healthResult.message());
  }

  public UnauthenticatedApi getUnauthenticatedApi() {
    return new UnauthenticatedApi(workspaceClient.getUnauthorizedApiClient());
  }

  public ResourceApi getResourceApi(String accessToken) {
    return new ResourceApi(workspaceClient.getApiClient(accessToken));
  }

  public ControlledAzureResourceApi getControlledAzureResourceApi(String accessToken) {
    return new ControlledAzureResourceApi(workspaceClient.getApiClient(accessToken));
  }

  public UUID getWorkspaceStorageResourceId(UUID workspaceId, String accessToken) {
    ResourceList resourceList =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                getResourceApi(accessToken)
                    .enumerateResources(
                        workspaceId,
                        null,
                        null,
                        ResourceType.AZURE_STORAGE_CONTAINER,
                        StewardshipType.CONTROLLED));
    // assuming for now that there's only one AZURE_STORAGE_CONTAINER resource per workspace
    return resourceList.getResources().get(0).getMetadata().getResourceId();
  }

  public String writePermissionString = "w";

  public String getSasTokenForFile(
      UUID workspaceId, String blobName, String sasPermissions, String accessToken) {
    UUID resourceId = getWorkspaceStorageResourceId(workspaceId, accessToken);
    Long sasExpirationDuration = 24 * 60 * 60L; // 24 hours in seconds; 24h is the max allowed
    // createAzureStorageContainerSasToken(UUID workspaceId, UUID resourceId, String sasIpRange,
    // Long sasExpirationDuration, String sasPermissions, String sasBlobName)
    CreatedAzureStorageContainerSasToken createdAzureStorageContainerSasToken =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                getControlledAzureResourceApi(accessToken)
                    .createAzureStorageContainerSasToken(
                        workspaceId,
                        resourceId,
                        null, // sasIpRange ???
                        sasExpirationDuration,
                        sasPermissions,
                        blobName));
    return createdAzureStorageContainerSasToken.getUrl();
  }

  interface WorkspaceAction<T> {
    T execute() throws ApiException;
  }

  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, WorkspaceAction<T> action)
      throws WorkspaceServiceApiException {

    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new WorkspaceServiceApiException(e);
          }
        });
  }
}
