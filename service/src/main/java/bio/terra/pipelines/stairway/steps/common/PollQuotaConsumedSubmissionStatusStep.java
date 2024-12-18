package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.app.configuration.internal.WdlPipelineConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.steps.utils.RawlsSubmissionStepHelper;
import bio.terra.stairway.*;
import java.util.UUID;
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
  private final WdlPipelineConfiguration wdlPipelineConfiguration;
  private final Logger logger =
      LoggerFactory.getLogger(PollQuotaConsumedSubmissionStatusStep.class);

  public PollQuotaConsumedSubmissionStatusStep(
      RawlsService rawlsService,
      SamService samService,
      WdlPipelineConfiguration wdlPipelineConfiguration) {
    this.samService = samService;
    this.rawlsService = rawlsService;
    this.wdlPipelineConfiguration = wdlPipelineConfiguration;
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

    RawlsSubmissionStepHelper rawlsSubmissionStepHelper =
        new RawlsSubmissionStepHelper(
            rawlsService, samService, controlWorkspaceProject, controlWorkspaceName, logger);
    return rawlsSubmissionStepHelper.pollRawlsSubmissionHelper(
        quotaSubmissionId, wdlPipelineConfiguration.getQuotaConsumedPollingIntervalSeconds());
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about polling a cromwell submission
    return StepResult.getStepResultSuccess();
  }
}
