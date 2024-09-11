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
  private static final String DATA_TABLE_REFERENCE_PREFIX = "this.";

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
   * validates a method config against the expected version, data table entity, inputs, and outputs.
   * It does this by checking that the wdl method version matches what we expect and that all expected
   * inputs/outputs both exist in the method config's list of inputs/outputs and that the data table
   * reference for each input/output are what we expect.
   *
   * @param methodConfiguration - method config to validate against
   * @param dataTableEntityName - data table entity name to use, should be the pipeline name
   * @param wdlWorkflowName - name of the wdl workflow, used to construct the full wdl variable name
   * @param pipelineInputDefinitions - list of input definitions for a pipeline
   * @param pipelineOutputDefinitions - list of output definitions for a pipeline
   * @param wdlMethodVersion - version of wdl we should be submitting
   * @return - whether the method config matches what we expect with a boolean
   */
  public boolean validateMethodConfig(
      MethodConfiguration methodConfiguration,
      String dataTableEntityName,
      String wdlWorkflowName,
      List<PipelineInputDefinition> pipelineInputDefinitions,
      List<PipelineOutputDefinition> pipelineOutputDefinitions,
      String wdlMethodVersion) {
    // validate wdl method version
    if (!methodConfiguration.getMethodRepoMethod().getMethodVersion().equals(wdlMethodVersion)) {
      return false;
    }
    // validate data table entity name
    if (!methodConfiguration.getRootEntityType().equals(dataTableEntityName)) {
      return false;
    }
    // validate inputs
    HashMap<?, ?> methodConfigInputs =
        objectMapper.convertValue(methodConfiguration.getInputs(), HashMap.class);
    for (PipelineInputDefinition pipelineInputDefinition : pipelineInputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineInputDefinition.getWdlVariableName();
      String fullDataTableReference =
          DATA_TABLE_REFERENCE_PREFIX + pipelineInputDefinition.getWdlVariableName();
      if (!methodConfigInputs.containsKey(fullWdlVariableName)
          || !methodConfigInputs.get(fullWdlVariableName).equals(fullDataTableReference)) {
        return false;
      }
    }
    // validate outputs
    HashMap<?, ?> methodConfigOutputs =
        objectMapper.convertValue(methodConfiguration.getOutputs(), HashMap.class);
    for (PipelineOutputDefinition pipelineOutputDefinition : pipelineOutputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineOutputDefinition.getWdlVariableName();
      String fullDataTableReference =
          DATA_TABLE_REFERENCE_PREFIX + pipelineOutputDefinition.getWdlVariableName();
      if (!methodConfigOutputs.containsKey(fullWdlVariableName)
          || !methodConfigOutputs.get(fullWdlVariableName).equals(fullDataTableReference)) {
        return false;
      }
    }

    // all validation checks passed, return true
    return true;
  }

  /**
   * this method takes a pre-existing method configuration and updates it to match what the service
   * expects. It does this by setting the input and outputs fields to match what we expect both from
   * what keys we expect and the data table reference we expect.
   *
   * <p>It also sets the wdl method version to the version we should be running as well as updating
   * the methodRepoMethod method uri by swapping the old version string with the expected version
   * string.
   *
   * @param methodConfiguration - method configuration to update
   * @param dataTableEntityName - data table entity name to use, should be the pipeline name
   * @param wdlWorkflowName - name of the wdl workflow, used to construct the full wdl variable name
   * @param pipelineInputDefinitions - list of input definitions for a pipeline
   * @param pipelineOutputDefinitions - list of output definitions for a pipeline
   * @param wdlMethodVersion - version of wdl we should be submitting
   * @return - a method configuration that has its methodRepoMethod's methodVersion update and its
   *     inputs and outputs
   */
  public MethodConfiguration updateMethodConfigToBeValid(
      MethodConfiguration methodConfiguration,
      String dataTableEntityName,
      String wdlWorkflowName,
      List<PipelineInputDefinition> pipelineInputDefinitions,
      List<PipelineOutputDefinition> pipelineOutputDefinitions,
      String wdlMethodVersion) {
    // set wdl method version and uri to match the expected version
    String oldVersion = methodConfiguration.getMethodRepoMethod().getMethodVersion();
    String oldUri = methodConfiguration.getMethodRepoMethod().getMethodUri();
    methodConfiguration.getMethodRepoMethod().setMethodVersion(wdlMethodVersion);
    methodConfiguration
        .getMethodRepoMethod()
        .setMethodUri(oldUri.replace(oldVersion, wdlMethodVersion));

    // set data table entity
    methodConfiguration.setRootEntityType(dataTableEntityName);

    // set inputs
    HashMap<String, String> expectedInputs = new HashMap<>();
    for (PipelineInputDefinition pipelineInputDefinition : pipelineInputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineInputDefinition.getWdlVariableName();
      String fullDataTableReference =
          DATA_TABLE_REFERENCE_PREFIX + pipelineInputDefinition.getWdlVariableName();
      expectedInputs.put(fullWdlVariableName, fullDataTableReference);
    }
    methodConfiguration.setInputs(expectedInputs);

    // set outputs
    HashMap<String, String> expectedOutputs = new HashMap<>();
    for (PipelineOutputDefinition pipelineOutputDefinition : pipelineOutputDefinitions) {
      String fullWdlVariableName =
          wdlWorkflowName + "." + pipelineOutputDefinition.getWdlVariableName();
      String fullDataTableReference =
          DATA_TABLE_REFERENCE_PREFIX + pipelineOutputDefinition.getWdlVariableName();
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
