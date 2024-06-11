package bio.terra.pipelines.dependencies.workspacemanager;

import bio.terra.workspace.client.ApiException;

public class WorkspaceServiceApiException extends WorkspaceServiceException {

  public WorkspaceServiceApiException(ApiException exception) {
    super("Workspace Manager returned an unsuccessful status code", exception);
  }

  public WorkspaceServiceApiException(String message) {
    super("Workspace Manager returned an unsuccessful status code");
  }
}
