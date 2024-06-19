package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.workspace.client.ApiException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class WorkspaceManagerServiceApiException extends WorkspaceManagerServiceException {

  public WorkspaceManagerServiceApiException(ApiException exception) {
    super("Workspace Manager returned an unsuccessful status code", exception);
  }

  public WorkspaceManagerServiceApiException(String message) {
    super("Workspace Manager returned an unsuccessful status code");
  }
}
