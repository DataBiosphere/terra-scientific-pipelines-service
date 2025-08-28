package bio.terra.pipelines.stairway.steps.common;

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
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.rawls.model.SubmissionReport;
import bio.terra.rawls.model.SubmissionRequest;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import java.math.BigDecimal;
import java.util.List;
import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step submits a wdl to cromwell using the rawls submission endpoint. The wdl that is run
 * depends on the workspace name and billing project provided to the step.
 *
 * <p>this step expects nothing from the working map
 *
 * <p>this step writes the submission_id to the working map using the submissionIdKey
 */
public class SubmitCromwellSubmissionStep implements Step {
  private final SamService samService;
  private final RawlsService rawlsService;

  private final String toolConfigKey;
  private final String submissionIdKey;

  private final Logger logger = LoggerFactory.getLogger(SubmitCromwellSubmissionStep.class);

  public SubmitCromwellSubmissionStep(
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
        ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT,
        ImputationJobMapKeys.CONTROL_WORKSPACE_NAME,
        toolConfigKey);

    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
    String controlWorkspaceName =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_NAME, String.class);
    String controlWorkspaceProject =
        inputParameters.get(ImputationJobMapKeys.CONTROL_WORKSPACE_BILLING_PROJECT, String.class);

    String description = inputParameters.get(JobMapKeys.DESCRIPTION, String.class);
    ToolConfig toolConfig = inputParameters.get(toolConfigKey, new TypeReference<>() {});
    logger.info(
        "Submitting method: {}, version: {}", toolConfig.methodName(), toolConfig.methodVersion());

    String methodName = toolConfig.methodName();
    String methodVersion = toolConfig.methodVersion();
    BigDecimal memoryRetryMultiplier = toolConfig.memoryRetryMultiplier();
    List<PipelineInputDefinition> inputDefinitions = toolConfig.inputDefinitions();
    List<PipelineOutputDefinition> outputDefinitions = toolConfig.outputDefinitions();
    boolean useCallCache = toolConfig.callCache();
    boolean deleteIntermediateOutputFiles = toolConfig.deleteIntermediateOutputFiles();
    boolean useReferenceDisks = toolConfig.useReferenceDisks();

    RawlsSubmissionStepHelper rawlsSubmissionStepHelper =
        new RawlsSubmissionStepHelper(
            rawlsService, samService, controlWorkspaceProject, controlWorkspaceName, logger);

    Optional<StepResult> validationResponse =
        rawlsSubmissionStepHelper.validateRawlsSubmissionMethodHelper(
            methodName, methodVersion, inputDefinitions, outputDefinitions, pipelineName);

    // if there is a validation response that means the validation failed so return it
    if (validationResponse.isPresent()) {
      return validationResponse.get();
    }

    // create submission request
    SubmissionRequest submissionRequest =
        new SubmissionRequest()
            .entityName(flightContext.getFlightId())
            .entityType(pipelineName.getValue())
            .useCallCache(useCallCache)
            .deleteIntermediateOutputFiles(deleteIntermediateOutputFiles)
            .useReferenceDisks(useReferenceDisks)
            .memoryRetryMultiplier(memoryRetryMultiplier)
            .userComment(
                "%s - %s for flight id: %s; description: %s"
                    .formatted(pipelineName, methodName, flightContext.getFlightId(), description))
            .methodConfigurationNamespace(controlWorkspaceProject)
            .methodConfigurationName(methodName);

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
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(submissionIdKey, submissionReport.getSubmissionId());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about submitting a run set
    return StepResult.getStepResultSuccess();
  }
}
