package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.pipelines.app.configuration.external.WorkspaceManagerServerConfiguration;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
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

  public ApiSystemStatusSystems checkHealthApiSystemStatus() {
    Result healthResult = checkHealth();
    return new ApiSystemStatusSystems()
        .ok(healthResult.isOk())
        .addMessagesItem(healthResult.message());
  }

  public UUID getWorkspaceStorageResourceId(UUID workspaceId, String accessToken) {
    ResourceList resourceList =
        executionWithRetryTemplate(
            listenerResetRetryTemplate,
            () ->
                workspaceManagerClient
                    .getResourceApi(accessToken)
                    .enumerateResources(
                        workspaceId,
                        null,
                        null,
                        ResourceType.AZURE_STORAGE_CONTAINER,
                        StewardshipType.CONTROLLED));

    // there should be only one AZURE_STORAGE_CONTAINER resource per workspace
    if (resourceList.getResources().size() != 1) {
      logger.warn(
          "Workspace {} has {} AZURE_STORAGE_CONTAINER resources, expected 1",
          workspaceId,
          resourceList.getResources().size());
    }

    return resourceList.getResources().get(0).getMetadata().getResourceId();
  }

  public static final String READ_PERMISSION_STRING = "r";

  public String getSasTokenForFile(
      UUID workspaceId, String fullFilePath, String sasPermissions, String accessToken) {
    logger.info("Calling WSM to get storage container id for workspace: {}", workspaceId);
    UUID resourceId = getWorkspaceStorageResourceId(workspaceId, accessToken);
    Long sasExpirationDuration =
        workspaceManagerServerConfiguration.sasExpirationDurationHours() * 60 * 60;

    // extract the blob name from the full file path
    String blobName = getBlobNameFromHttpUrl(fullFilePath, workspaceId);

    logger.info(
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

  protected String getBlobNameFromHttpUrl(String blobHttpUrl, UUID workspaceId) {
    // return the remaining string after the workspaceId
    return blobHttpUrl.substring(
        blobHttpUrl.indexOf(workspaceId.toString()) + workspaceId.toString().length() + 1);
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
}
