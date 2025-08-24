package bio.terra.pipelines.stairway.steps.imputation.gcp;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.pipelines.stairway.steps.utils.RawlsSubmissionStepHelper;
import bio.terra.rawls.model.SubmissionReport;
import bio.terra.rawls.model.SubmissionRequest;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.List;
import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step submits a submission to cromwell using the rawls submission endpoint. It uses
 * ImputationConfiguration in order to set some config options on the cromwell submission
 *
 * <p>this step expects nothing from the working map
 *
 * <p>this step writes submission_id to the working map
 */
public class SubmitImputationSubmissionStep implements Step {
  private final SamService samService;
  private final RawlsService rawlsService;
  private final ImputationConfiguration imputationConfiguration;

  private final Logger logger = LoggerFactory.getLogger(SubmitImputationSubmissionStep.class);

  public SubmitImputationSubmissionStep(
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
        JobMapKeys.PIPELINE_NAME,
        JobMapKeys.DESCRIPTION,
        ImputationJobMapKeys.WDL_METHOD_NAME,
        ImputationJobMapKeys.WDL_METHOD_VERSION,
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
        ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS);

    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);
    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
    String description = inputParameters.get(JobMapKeys.DESCRIPTION, String.class);
    String wdlMethodName = inputParameters.get(ImputationJobMapKeys.WDL_METHOD_NAME, String.class);
    String wdlMethodVersion =
        inputParameters.get(ImputationJobMapKeys.WDL_METHOD_VERSION, String.class);

    List<PipelineInputDefinition> inputDefinitions =
        inputParameters.get(
            ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS, new TypeReference<>() {});
    List<PipelineOutputDefinition> outputDefinitions =
        inputParameters.get(
            ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS, new TypeReference<>() {});

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();

    RawlsSubmissionStepHelper rawlsSubmissionStepHelper =
        new RawlsSubmissionStepHelper(
            rawlsService, samService, controlWorkspaceProject, controlWorkspaceName, logger);

    Optional<StepResult> validationResponse =
        rawlsSubmissionStepHelper.validateRawlsSubmissionMethodHelper(
            wdlMethodName, wdlMethodVersion, inputDefinitions, outputDefinitions, pipelineName);

    // if there is a validation response that means the validation failed so return it
    if (validationResponse.isPresent()) {
      return validationResponse.get();
    }

    // create submission request
    SubmissionRequest submissionRequest =
        new SubmissionRequest()
            .entityName(flightContext.getFlightId())
            .entityType(pipelineName.getValue())
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
    workingMap.put(ImputationJobMapKeys.SUBMISSION_ID, submissionReport.getSubmissionId());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about submitting a run set
    return StepResult.getStepResultSuccess();
  }
}
