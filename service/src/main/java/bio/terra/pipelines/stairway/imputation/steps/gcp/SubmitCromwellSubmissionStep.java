package bio.terra.pipelines.stairway.imputation.steps.gcp;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.rawls.model.SubmissionReport;
import bio.terra.rawls.model.SubmissionRequest;
import bio.terra.stairway.*;

/**
 * This step submits a submission to cromwell using the rawls submission endpoint. It uses
 * ImputationConfiguration in order to set some config options on the cromwell submission
 *
 * <p>this step expects nothing from the working map
 *
 * <p>this step writes submission_id to the working map
 */
public class SubmitCromwellSubmissionStep implements Step {
  private final SamService samService;
  private final RawlsService rawlsService;
  private final ImputationConfiguration imputationConfiguration;

  public SubmitCromwellSubmissionStep(
      RawlsService rawlsService,
      SamService samService,
      ImputationConfiguration imputationConfiguration) {
    this.samService = samService;
    this.rawlsService = rawlsService;
    this.imputationConfiguration = imputationConfiguration;
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
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        JobMapKeys.DESCRIPTION.getKeyName(),
        RunImputationJobFlightMapKeys.WDL_METHOD_NAME,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME);

    String controlWorkspaceName =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(
            RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);
    String description = inputParameters.get(JobMapKeys.DESCRIPTION.getKeyName(), String.class);
    String wdlMethodName =
        inputParameters.get(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, String.class);

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
            .useCallCache(imputationConfiguration.isUseCallCaching())
            .deleteIntermediateOutputFiles(imputationConfiguration.isDeleteIntermediateFiles())
            .useReferenceDisks(imputationConfiguration.isUseReferenceDisk())
            .userComment(
                "%s (%s) - flight id: %s; description: %s"
                    .formatted(
                        pipelineName, wdlMethodName, flightContext.getFlightId(), description))
            .methodConfigurationNamespace(controlWorkspaceProject)
            .methodConfigurationName(wdlMethodName);

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
    workingMap.put(RunImputationJobFlightMapKeys.SUBMISSION_ID, submissionReport.getSubmissionId());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about submitting a run set
    return StepResult.getStepResultSuccess();
  }
}