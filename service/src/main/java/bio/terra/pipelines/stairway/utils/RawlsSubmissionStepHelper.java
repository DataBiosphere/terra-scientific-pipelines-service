package bio.terra.pipelines.stairway.utils;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.rawls.model.MethodConfiguration;
import bio.terra.rawls.model.Submission;
import bio.terra.rawls.model.Workflow;
import bio.terra.rawls.model.WorkflowStatus;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;

public class RawlsSubmissionStepHelper {

  private final SamService samService;
  private final RawlsService rawlsService;
  private final String controlWorkspaceProject;
  private final String controlWorkspaceName;
  private final Logger logger;

  public RawlsSubmissionStepHelper(
      RawlsService rawlsService,
      SamService samService,
      String controlWorkspaceProject,
      String controlWorkspaceName,
      Logger logger) {
    this.rawlsService = rawlsService;
    this.samService = samService;
    this.controlWorkspaceProject = controlWorkspaceProject;
    this.controlWorkspaceName = controlWorkspaceName;
    this.logger = logger;
  }

  public StepResult pollRawlsSubmissionHelper(UUID submissionId, Long secondsToSleep)
      throws InterruptedException {
    // poll until all runs are in a finalized state
    Submission submissionResponse = null;
    boolean stillRunning = true;
    try {
      while (stillRunning) {
        submissionResponse =
            rawlsService.getSubmissionStatus(
                samService.getTeaspoonsServiceAccountToken(),
                controlWorkspaceProject,
                controlWorkspaceName,
                submissionId);
        stillRunning = RawlsService.submissionIsRunning(submissionResponse);
        if (stillRunning) {
          logger.info("Polling Started, sleeping for {} seconds", secondsToSleep);
          TimeUnit.SECONDS.sleep(secondsToSleep);
        }
      }
    } catch (RawlsServiceApiException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    // if there are any non-successful workflows, fatally fail the step
    List<Workflow> failedRunLogs =
        submissionResponse.getWorkflows().stream()
            .filter(workflow -> !workflow.getStatus().equals(WorkflowStatus.SUCCEEDED))
            .toList();
    if (failedRunLogs.isEmpty()) {
      return StepResult.getStepResultSuccess();
    } else {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException(
              "Not all runs succeeded for submission: " + submissionId));
    }
  }

  public Optional<StepResult> validateRawlsSubmissionMethodHelper(
      String wdlMethodName,
      String wdlMethodVersion,
      List<PipelineInputDefinition> inputDefinitions,
      List<PipelineOutputDefinition> outputDefinitions,
      PipelinesEnum pipelineName) {
    MethodConfiguration methodConfiguration;
    try {
      // grab current method config and validate it
      methodConfiguration =
          rawlsService.getCurrentMethodConfigForMethod(
              samService.getTeaspoonsServiceAccountToken(),
              controlWorkspaceProject,
              controlWorkspaceName,
              wdlMethodName);
    } catch (RawlsServiceApiException e) {
      // if we fail to grab the method config then retry
      return Optional.of(new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e));
    }
    boolean validMethodConfig =
        rawlsService.validateMethodConfig(
            methodConfiguration,
            pipelineName.getValue(),
            wdlMethodName,
            inputDefinitions,
            outputDefinitions,
            wdlMethodVersion);

    // if not a valid method config, set the method config to what we think it should be.  This
    // shouldn't happen
    if (!validMethodConfig) {
      logger.warn(
          "found method config that was not valid for billing project: {}, workspace: {}, method name: {}, methodConfigVersion: {}",
          controlWorkspaceProject,
          controlWorkspaceName,
          wdlMethodName,
          methodConfiguration.getMethodConfigVersion());

      MethodConfiguration updatedMethodConfiguration =
          rawlsService.updateMethodConfigToBeValid(
              methodConfiguration,
              pipelineName.getValue(),
              wdlMethodName,
              inputDefinitions,
              outputDefinitions,
              wdlMethodVersion);
      try {
        // update method config version, inputs, and outputs
        rawlsService.setMethodConfigForMethod(
            samService.getTeaspoonsServiceAccountToken(),
            updatedMethodConfiguration,
            controlWorkspaceProject,
            controlWorkspaceName,
            wdlMethodName);
      } catch (RawlsServiceApiException e) {
        // if we fail to update the method config then retry
        return Optional.of(new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e));
      }
    }
    return Optional.empty();
  }
}
