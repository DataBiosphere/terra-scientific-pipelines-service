package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.imputation.ImputationJobMapKeys;
import bio.terra.rawls.model.SubmissionReport;
import bio.terra.rawls.model.SubmissionRequest;
import bio.terra.stairway.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step submits a submission of a quota wdl to cromwell using the rawls submission endpoint.
 * The quota consumed wdl that is run depends on the workspace name provided to the step.
 *
 * <p>this step expects nothing from the working map
 *
 * <p>this step writes quota_submission_id to the working map
 */
public class SubmitQuotaConsumedSubmissionStep implements Step {
  private final SamService samService;
  private final RawlsService rawlsService;

  private final Logger logger = LoggerFactory.getLogger(SubmitQuotaConsumedSubmissionStep.class);

  public SubmitQuotaConsumedSubmissionStep(RawlsService rawlsService, SamService samService) {
    this.samService = samService;
    this.rawlsService = rawlsService;
  }

  @Override
  @SuppressWarnings(
      "java:S2259") // suppress warning for possible NPE when calling pipelineName.getValue(),
  //  since we do validate that pipelineName is not null in `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME);

    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();

    // create submission request
    SubmissionRequest submissionRequest =
        new SubmissionRequest()
            .entityName(flightContext.getFlightId())
            .entityType(
                pipelineName.getValue()) // this must match the configuration the method is set to
            // launch with.  Will be addressed in
            // https://broadworkbench.atlassian.net/browse/TSPS-301
            .useCallCache(true)
            .deleteIntermediateOutputFiles(true)
            .useReferenceDisks(false)
            .userComment(
                "%s - getting quota consumed for flight id: %s"
                    .formatted(pipelineName, flightContext.getFlightId()))
            .methodConfigurationNamespace(controlWorkspaceProject)
            .methodConfigurationName("QuotaConsumed");

    // submit workflow to rawls
    SubmissionReport submissionReport;
    try {
      submissionReport =
          rawlsService.submitWorkflow(
              samService.getTeaspoonsServiceAccountToken(),
              submissionRequest,
              controlWorkspaceProject,
              controlWorkspaceName);
    } catch (RawlsServiceApiException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    // add submission id to working map to be used for polling in downstream step
    workingMap.put(ImputationJobMapKeys.QUOTA_SUBMISSION_ID, submissionReport.getSubmissionId());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about submitting a run set
    return StepResult.getStepResultSuccess();
  }
}
