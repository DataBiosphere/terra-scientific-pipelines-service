package bio.terra.pipelines.stairway;

import bio.terra.common.exception.InternalServerErrorException;
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
 * not in a final state, this will step will poll again after an interval of time. Once the
 * submission is finalized then it will see if the workflows are all successful and if so will
 * succeed otherwise will fail.
 *
 * <p>this step expects quota submission id to be provided in the working map
 */
public class PollQuotaConsumedSubmissionStatusStep implements Step {
  private final RawlsService rawlsService;
  private final SamService samService;
  private final Logger logger =
      LoggerFactory.getLogger(PollQuotaConsumedSubmissionStatusStep.class);

  public PollQuotaConsumedSubmissionStatusStep(RawlsService rawlsService, SamService samService) {
    this.samService = samService;
    this.rawlsService = rawlsService;
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
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.QUOTA_SUBMISSION_ID);

    UUID quotaSubmissionId = workingMap.get(ImputationJobMapKeys.QUOTA_SUBMISSION_ID, UUID.class);

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
                quotaSubmissionId);
        stillRunning = RawlsService.submissionIsRunning(submissionResponse);
        if (stillRunning) {
          logger.info("Polling Started, sleeping for {} seconds", 60);
          TimeUnit.SECONDS.sleep(60);
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
              "Not all runs succeeded for submission: " + quotaSubmissionId));
    }
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about polling a cromwell submission
    return StepResult.getStepResultSuccess();
  }
}
