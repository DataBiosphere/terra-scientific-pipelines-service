package bio.terra.pipelines.dependencies.rawls;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.rawls.client.ApiException;
import bio.terra.rawls.model.*;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.HashMap;
import java.util.List;
import java.util.UUID;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

@Service
public class RawlsService implements HealthCheck {
  private final RawlsClient rawlsClient;
  private final RetryTemplate listenerResetRetryTemplate;
  private final ObjectMapper objectMapper;

  private static final List<SubmissionStatus> FINAL_RUN_STATES =
      List.of(SubmissionStatus.ABORTED, SubmissionStatus.DONE);

  @Autowired
  public RawlsService(
      RawlsClient rawlsClient,
      RetryTemplate listenerResetRetryTemplate,
      ObjectMapper objectMapper) {
    this.rawlsClient = rawlsClient;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
    this.objectMapper = objectMapper;
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

  public WorkspaceDetails getWorkspaceDetails(
      String accessToken, String workspaceNamespace, String workspaceName) {
    List<String> fields = List.of("workspace.bucketName", "workspace.googleProject");
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getWorkspacesApi(accessToken)
                .listWorkspaceDetails(workspaceNamespace, workspaceName, fields)
                .getWorkspace());
  }

  public String getWorkspaceBucketName(WorkspaceDetails workspaceDetails) {
    String bucketName = workspaceDetails.getBucketName();
    if (bucketName == null) {
      throw new InternalServerErrorException("Workspace bucket name is not defined");
    }
    return bucketName;
  }

  public String getWorkspaceGoogleProject(WorkspaceDetails workspaceDetails) {
    String googleProject = workspaceDetails.getGoogleProject();
    if (googleProject == null) {
      throw new InternalServerErrorException("Workspace google project is not defined");
    }
    return googleProject;
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

  public Entity getDataTableEntity(
      String accessToken,
      String workspaceNamespace,
      String workspaceName,
      String entityType,
      String entityName) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getEntitiesApi(accessToken)
                .getEntity(workspaceNamespace, workspaceName, entityType, entityName, null, null));
  }

  public MethodConfiguration getCurrentMethodConfigForMethod(
      String accessToken, String workspaceNamespace, String workspaceName, String methodName) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getMethodConfigsApi(accessToken)
                .getMethodConfiguration(
                    workspaceNamespace, workspaceName, workspaceNamespace, methodName));
  }

  public ValidatedMethodConfiguration setMethodConfigForMethod(
      String accessToken,
      MethodConfiguration methodConfiguration,
      String workspaceNamespace,
      String workspaceName,
      String methodName) {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            rawlsClient
                .getMethodConfigsApi(accessToken)
                .updateMethodConfiguration(
                    methodConfiguration,
                    workspaceNamespace,
                    workspaceName,
                    workspaceNamespace,
                    methodName));
  }

  // returns true if submission is in a running state
  public static boolean submissionIsRunning(Submission submission) {
    return !FINAL_RUN_STATES.contains(submission.getStatus());
  }

  /**
   * validates a method config against the expected version, inputs, and outputs
   *
   * @param methodConfiguration - method config to validate against
   * @param wdlWorkflowName - name of the wdl workflow, used to construct the full wdl variable name
   * @param pipelineInputDefinitions - list of input definitions for a pipeline
   * @param pipelineOutputDefinitions - list of output definitions for a pipeline
   * @param wdlMethodVersion - version of wdl we should be submitting
   * @return - whether the method config matches what we expect with a boolean
   */
  public boolean validateMethodConfig(
      MethodConfiguration methodConfiguration,
      String wdlWorkflowName,
      List<PipelineInputDefinition> pipelineInputDefinitions,
      List<PipelineOutputDefinition> pipelineOutputDefinitions,
      String wdlMethodVersion) {
    boolean isValid = true;
    // check wdl method version
    if (!methodConfiguration.getMethodRepoMethod().getMethodVersion().equals(wdlMethodVersion)) {
      isValid = false;
    }
    // validate inputs
    HashMap<?, ?> methodConfigInputs =
        objectMapper.convertValue(methodConfiguration.getInputs(), HashMap.class);
    for (PipelineInputDefinition pipelineInputDefinition : pipelineInputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineInputDefinition.getWdlVariableName();
      String fullDataTableReference = "this." + pipelineInputDefinition.getWdlVariableName();
      if (!methodConfigInputs.containsKey(fullWdlVariableName)
          || !methodConfigInputs.get(fullWdlVariableName).equals(fullDataTableReference)) {
        isValid = false;
      }
    }
    // validate outputs
    HashMap<?, ?> methodConfigOutputs =
        objectMapper.convertValue(methodConfiguration.getOutputs(), HashMap.class);
    for (PipelineOutputDefinition pipelineOutputDefinition : pipelineOutputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineOutputDefinition.getWdlVariableName();
      String fullDataTableReference = "this." + pipelineOutputDefinition.getWdlVariableName();
      if (!methodConfigOutputs.containsKey(fullWdlVariableName)
          || !methodConfigOutputs.get(fullWdlVariableName).equals(fullDataTableReference)) {
        isValid = false;
      }
    }

    return isValid;
  }

  /**
   * this method takes a pre-existing method configuration and updates it to match what the service
   * expects. It does this by updating the methodRepoMethod methodVersion and by setting the inputs
   * and outputs fields to match the references we expect.
   *
   * @param methodConfiguration - method configuration to update
   * @param methodConfiguration - method config to validate against
   * @param wdlWorkflowName - name of the wdl workflow, used to construct the full wdl variable name
   * @param pipelineInputDefinitions - list of input definitions for a pipeline
   * @param pipelineOutputDefinitions - list of output definitions for a pipeline
   * @param wdlMethodVersion - version of wdl we should be submitting
   * @return - a method configuration that has its methodRepoMethod's methodVersion update and its
   *     inputs and outputs
   */
  public MethodConfiguration updateMethodConfigToBeValid(
      MethodConfiguration methodConfiguration,
      String wdlWorkflowName,
      List<PipelineInputDefinition> pipelineInputDefinitions,
      List<PipelineOutputDefinition> pipelineOutputDefinitions,
      String wdlMethodVersion) {
    // update wdl method version and uri to match version
    String oldVersion = methodConfiguration.getMethodRepoMethod().getMethodVersion();
    String oldUri = methodConfiguration.getMethodRepoMethod().getMethodUri();
    methodConfiguration.getMethodRepoMethod().setMethodVersion(wdlMethodVersion);
    methodConfiguration
        .getMethodRepoMethod()
        .setMethodUri(oldUri.replace(oldVersion, wdlMethodVersion));

    // update inputs
    HashMap<String, String> expectedInputs = new HashMap<>();
    for (PipelineInputDefinition pipelineInputDefinition : pipelineInputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineInputDefinition.getWdlVariableName();
      String fullDataTableReference = "this." + pipelineInputDefinition.getWdlVariableName();
      expectedInputs.put(fullWdlVariableName, fullDataTableReference);
    }
    methodConfiguration.setInputs(expectedInputs);

    // validate outputs
    HashMap<String, String> expectedOutputs = new HashMap<>();
    for (PipelineOutputDefinition pipelineOutputDefinition : pipelineOutputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineOutputDefinition.getWdlVariableName();
      String fullDataTableReference = "this." + pipelineOutputDefinition.getWdlVariableName();
      expectedOutputs.put(fullWdlVariableName, fullDataTableReference);
    }
    methodConfiguration.setOutputs(expectedOutputs);
    return methodConfiguration;
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
