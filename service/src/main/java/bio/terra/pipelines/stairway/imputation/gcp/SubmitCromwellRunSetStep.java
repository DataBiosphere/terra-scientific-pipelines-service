package bio.terra.pipelines.stairway.imputation.gcp;

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
 * This step submits a run set to cromwell using the rawls submission endpoint.
 *
 * <p>this step expects nothing from the working map
 *
 * <p>this step writes submission_id to the working map
 */
public class SubmitCromwellRunSetStep implements Step {
  private final SamService samService;
  private final RawlsService rawlsService;

  public SubmitCromwellRunSetStep(RawlsService rawlsService, SamService samService) {
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
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_PROJECT,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME);

    String controlWorkspaceName =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_PROJECT, String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();

    // create submission request
    SubmissionRequest submissionRequest =
        new SubmissionRequest()
            .entityName(flightContext.getFlightId())
            .entityType(
                pipelineName.getValue()) // this must match the configuration the method is set to
            // launch with.  Do we want to modify the current method
            // configuration each time so this matches?
            .useCallCache(true)
            .deleteIntermediateOutputFiles(false)
            .useReferenceDisks(false)
            .userComment("this is my flightid: %s".formatted(flightContext.getFlightId()))
            .methodConfigurationNamespace(controlWorkspaceProject)
            .methodConfigurationName(
                "ImputationBeagleEmpty"); // need to not hardcode this eventually

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
