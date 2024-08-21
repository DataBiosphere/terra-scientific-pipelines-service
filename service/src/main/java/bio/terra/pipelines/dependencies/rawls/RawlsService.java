package bio.terra.pipelines.dependencies.rawls;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.rawls.client.ApiException;
import bio.terra.rawls.model.*;
import java.util.List;
import java.util.UUID;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class RawlsService implements HealthCheck {
  private final RawlsClient rawlsClient;
  private final RetryTemplate listenerResetRetryTemplate;

  private static final List<SubmissionStatus> FINAL_RUN_STATES =
      List.of(SubmissionStatus.ABORTED, SubmissionStatus.DONE);

  public RawlsService(RawlsClient rawlsClient, RetryTemplate listenerResetRetryTemplate) {
    this.rawlsClient = rawlsClient;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  @Override
  public Result checkHealth() {
    // systemStatus is a void method, throws ApiException if status is not healthy
    try {
      rawlsClient.getStatusApi().systemStatus();
      return new Result(true, "Rawls is ok");
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  public SubmissionReport submitWorkflow(
      String accessToken,
      SubmissionRequest submissionRequest,
      String workspaceNamespace,
      String workspaceName) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getSubmissionsApi(accessToken)
                .createSubmission(submissionRequest, workspaceNamespace, workspaceName));
  }

  public Submission getSubmissionStatus(
      String accessToken, String workspaceNamespace, String workspaceName, UUID submissionId) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getSubmissionsApi(accessToken)
                .getSubmissionStatus(workspaceNamespace, workspaceName, submissionId.toString()));
  }

  public Entity upsertDataTableEntity(
      String accessToken, String workspaceNamespace, String workspaceName, Entity entity) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getEntitiesApi(accessToken)
                .createEntity(entity, workspaceNamespace, workspaceName));
  }

  // returns true if submission is in a running state
  public static boolean submissionIsRunning(Submission submission) {
    return !FINAL_RUN_STATES.contains(submission.getStatus());
  }

  interface RawlsAction<T> {
    T execute() throws ApiException;
  }

  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, RawlsAction<T> action)
      throws RawlsServiceApiException {

    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new RawlsServiceApiException(e);
          }
        });
  }
}
