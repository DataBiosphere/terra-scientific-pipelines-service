package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.steps.utils.RawlsSubmissionStepHelper;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step polls rawls for a submission until all runs are in a finalized state. If submission is
 * not in a final state, this will step will poll again after an interval of time. Once the
 * submission is finalized then it will see if the workflows are all successful and if so will
 * succeed otherwise will fail.
 *
 * <p>this step expects the submission id to be provided in the working map
 */
public class PollCromwellSubmissionStatusStep implements Step {
  private final RawlsService rawlsService;
  private final SamService samService;
  private final String toolConfigKey;
  private final String submissionIdKey;
  private final Logger logger = LoggerFactory.getLogger(PollCromwellSubmissionStatusStep.class);

  public PollCromwellSubmissionStatusStep(
      RawlsService rawlsService,
      SamService samService,
      String toolConfigKey,
      String submissionIdKey) {
    this.samService = samService;
    this.rawlsService = rawlsService;
    this.toolConfigKey = toolConfigKey;
    this.submissionIdKey = submissionIdKey;
  }

  @Override
  @SuppressWarnings("java:S2259") // suppress warning for possible NPE when calling
  // toolConfig.pollingIntervalSeconds(),
  //  since we do validate that toolConfig is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) throws InterruptedException {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        toolConfigKey);

    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    ToolConfig toolConfig = inputParameters.get(toolConfigKey, new TypeReference<>() {});

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, submissionIdKey);

    UUID quotaSubmissionId = workingMap.get(submissionIdKey, UUID.class);

    Long pollingIntervalSeconds = toolConfig.pollingIntervalSeconds();
    logger.info(
        "Polling Rawls for {} submission {} status", toolConfig.methodName(), quotaSubmissionId);
    RawlsSubmissionStepHelper rawlsSubmissionStepHelper =
        new RawlsSubmissionStepHelper(
            rawlsService, samService, controlWorkspaceProject, controlWorkspaceName, logger);
    return rawlsSubmissionStepHelper.pollRawlsSubmissionHelper(
        quotaSubmissionId, pollingIntervalSeconds);
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about polling a cromwell submission
    return StepResult.getStepResultSuccess();
  }
}
