package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.pipelines.app.configuration.external.WorkspaceManagerServerConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.workspace.client.ApiException;
import bio.terra.workspace.model.CreatedAzureStorageContainerSasToken;
import bio.terra.workspace.model.ResourceList;
import bio.terra.workspace.model.ResourceType;
import bio.terra.workspace.model.StewardshipType;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class WorkspaceManagerService implements HealthCheck {
  private final WorkspaceManagerClient workspaceManagerClient;
  private final WorkspaceManagerServerConfiguration workspaceManagerServerConfiguration;
  private final RetryTemplate listenerResetRetryTemplate;
  private static final Logger logger = LoggerFactory.getLogger(WorkspaceManagerService.class);

  public WorkspaceManagerService(
      WorkspaceManagerClient workspaceManagerClient,
      WorkspaceManagerServerConfiguration workspaceManagerServerConfiguration,
      RetryTemplate listenerResetRetryTemplate) {
    this.workspaceManagerClient = workspaceManagerClient;
    this.workspaceManagerServerConfiguration = workspaceManagerServerConfiguration;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  public Result checkHealth() {
    // No access token needed since this is an unauthenticated API.
    try {
      workspaceManagerClient
          .getUnauthenticatedApi()
          .serviceStatus(); // Workspace Manager's serviceStatus() is a void method
      return new Result(true, "Workspace Manager is ok");
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  private ResourceList getWorkspaceResources(
      UUID workspaceId,
      ResourceType resourceType,
      StewardshipType stewardshipType,
      String accessToken) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            workspaceManagerClient
                .getResourceApi(accessToken)
                .enumerateResources(workspaceId, null, null, resourceType, stewardshipType));
  }

  private String getSasUrlForBlob(
      UUID workspaceId,
      UUID resourceId,
      String blobName,
      Long sasExpirationDuration,
      String sasPermissions,
      String accessToken) {

    logger.debug(
        "Calling WSM to get SAS token (permissions: {}) for blob: {}", sasPermissions, blobName);

    CreatedAzureStorageContainerSasToken createdAzureStorageContainerSasToken =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                workspaceManagerClient
                    .getControlledAzureResourceApi(accessToken)
                    .createAzureStorageContainerSasToken(
                        workspaceId,
                        resourceId,
                        null,
                        sasExpirationDuration,
                        sasPermissions,
                        blobName));
    return createdAzureStorageContainerSasToken.getUrl();
  }

  interface WorkspaceAction<T> {
    T execute() throws ApiException;
  }

  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, WorkspaceAction<T> action)
      throws WorkspaceManagerServiceApiException {

    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new WorkspaceManagerServiceApiException(e);
          }
        });
  }

  public UUID getWorkspaceStorageResourceId(UUID workspaceId, String accessToken) {
    ResourceList resourceList =
        getWorkspaceResources(
            workspaceId,
            ResourceType.AZURE_STORAGE_CONTAINER,
            StewardshipType.CONTROLLED,
            accessToken);

    // there should be only one AZURE_STORAGE_CONTAINER resource per workspace
    if (resourceList.getResources().size() != 1) {
      logger.warn(
          "Workspace {} has {} AZURE_STORAGE_CONTAINER resources, expected 1",
          workspaceId,
          resourceList.getResources().size());
    }

    return resourceList.getResources().get(0).getMetadata().getResourceId();
  }

  public String getWriteSasUrlForBlob(
      UUID workspaceId, UUID resourceId, String blobName, String accessToken) {
    Long sasExpirationDurationSeconds =
        workspaceManagerServerConfiguration.sasExpirationDurationHours() * 60 * 60;
    return getSasUrlForBlob(
        workspaceId, resourceId, blobName, sasExpirationDurationSeconds, "w", accessToken);
  }

  public String getReadSasUrlForBlob(
      UUID workspaceId, UUID resourceId, String blobName, String accessToken) {
    Long sasExpirationDurationSeconds =
        workspaceManagerServerConfiguration.sasExpirationDurationHours() * 60 * 60;
    return getSasUrlForBlob(
        workspaceId, resourceId, blobName, sasExpirationDurationSeconds, "r", accessToken);
  }
}
