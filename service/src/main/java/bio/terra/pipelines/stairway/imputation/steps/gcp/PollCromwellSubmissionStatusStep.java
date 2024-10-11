package bio.terra.pipelines.stairway.imputation.steps.gcp;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.imputation.ImputationJobMapKeys;
import bio.terra.rawls.model.Submission;
import bio.terra.rawls.model.Workflow;
import bio.terra.rawls.model.WorkflowStatus;
import bio.terra.stairway.*;
import java.util.List;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step polls rawls for a submission until all runs are in a finalized state. If submission is
 * not in a final state, this will step will wait
 * bio.terra.pipelines.app.configuration.internal.ImputationConfiguration#cromwellSubmissionPollingIntervalInSeconds
 * seconds before polling again. Once the submission is finalized then it will see if the workflows
 * are all successful and if so will succeed otherwise will fail.
 *
 * <p>this step expects submission id to be provided in the working map
 */
public class PollCromwellSubmissionStatusStep implements Step {
  private final RawlsService rawlsService;
  private final SamService samService;
  private final ImputationConfiguration imputationConfiguration;
  private final Logger logger = LoggerFactory.getLogger(PollCromwellSubmissionStatusStep.class);

  public PollCromwellSubmissionStatusStep(
      SamService samService,
      RawlsService rawlsService,
      ImputationConfiguration imputationConfiguration) {
    this.samService = samService;
    this.rawlsService = rawlsService;
    this.imputationConfiguration = imputationConfiguration;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) throws InterruptedException {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT);
    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.SUBMISSION_ID);

    UUID submissionId = workingMap.get(ImputationJobMapKeys.SUBMISSION_ID, UUID.class);

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
          logger.info(
              "Polling Started, sleeping for {} seconds",
              imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
          TimeUnit.SECONDS.sleep(
              imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
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

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about polling a cromwell submission
    return StepResult.getStepResultSuccess();
  }
}
