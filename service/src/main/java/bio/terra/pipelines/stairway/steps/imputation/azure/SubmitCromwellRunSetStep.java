package bio.terra.pipelines.stairway.steps.imputation.azure;

import bio.terra.cbas.model.*;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.*;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.List;
import java.util.UUID;

/**
 * This step submits a run set to cromwell. It first finds the correct method version id to run
 * based on the method name passed to it. It then generates a run set request that links that
 * request to wds record(s). It then submits that request to cbas
 *
 * <p>this step expects the cbas uri to be passed in through the working map
 *
 * <p>this step writes run set id to the working map
 */
public class SubmitCromwellRunSetStep implements Step {
  private final CbasService cbasService;
  private final SamService samService;
  private final PipelinesService pipelinesService;
  private final CbasConfiguration cbasConfiguration;

  public SubmitCromwellRunSetStep(
      CbasService cbasService,
      SamService samService,
      PipelinesService pipelinesService,
      CbasConfiguration cbasConfiguration) {
    this.cbasService = cbasService;
    this.samService = samService;
    this.pipelinesService = pipelinesService;
    this.cbasConfiguration = cbasConfiguration;
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
        JobMapKeys.DESCRIPTION,
        JobMapKeys.PIPELINE_NAME,
        ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
        ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS,
        ImputationJobMapKeys.WDL_METHOD_NAME);

    String description = inputParameters.get(JobMapKeys.DESCRIPTION, String.class);
    PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
    List<PipelineInputDefinition> allInputDefinitions =
        inputParameters.get(
            ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS, new TypeReference<>() {});
    List<PipelineOutputDefinition> outputDefinitions =
        inputParameters.get(
            ImputationJobMapKeys.PIPELINE_OUTPUT_DEFINITIONS, new TypeReference<>() {});
    String wdlMethodName = inputParameters.get(ImputationJobMapKeys.WDL_METHOD_NAME, String.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.CBAS_URI);

    String cbasUri = workingMap.get(ImputationJobMapKeys.CBAS_URI, String.class);

    // grab methodVersionId needed to submit a submission
    MethodListResponse methodListResponse =
        cbasService.getAllMethods(cbasUri, samService.getTeaspoonsServiceAccountToken());
    UUID methodVersionId =
        CbasService.getMethodVersionIdFromMethodListResponse(methodListResponse, wdlMethodName);
    if (methodVersionId == null) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException(
              "couldn't find method version id for method name " + wdlMethodName));
    }

    // prepare cbas submission (run set) object
    // outputs are hardcoded for now (to be constructed dynamically in TSPS-197)
    RunSetRequest runSetRequest =
        new RunSetRequest()
            .runSetDescription(
                "%s (%s) - flight id: %s; description: %s"
                    .formatted(
                        pipelineName, wdlMethodName, flightContext.getFlightId(), description))
            .runSetName("%s - flightId %s".formatted(wdlMethodName, flightContext.getFlightId()))
            .methodVersionId(methodVersionId)
            .callCachingEnabled(cbasConfiguration.getCallCache())
            // define the WDS record to link to for inputs and outputs
            .wdsRecords(
                new WdsRecordSet()
                    .recordType(pipelineName.getValue())
                    .addRecordIdsItem(flightContext.getFlightId()));

    // add inputs
    List<WorkflowInputDefinition> cbasWorkflowInputDefinitions =
        pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(
            allInputDefinitions, wdlMethodName);
    for (WorkflowInputDefinition workflowInputDefinition : cbasWorkflowInputDefinitions) {
      runSetRequest.addWorkflowInputDefinitionsItem(workflowInputDefinition);
    }

    // add outputs
    List<WorkflowOutputDefinition> cbasWorkflowOutputDefinitions =
        pipelinesService.prepareCbasWorkflowOutputRecordUpdateDefinitions(
            outputDefinitions, wdlMethodName);
    for (WorkflowOutputDefinition workflowOutputDefinition : cbasWorkflowOutputDefinitions) {
      runSetRequest.addWorkflowOutputDefinitionsItem(workflowOutputDefinition);
    }

    // launch the submission
    RunSetStateResponse runSetStateResponse;
    try {
      runSetStateResponse =
          cbasService.createRunSet(
              cbasUri, samService.getTeaspoonsServiceAccountToken(), runSetRequest);
    } catch (CbasServiceApiException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }
    workingMap.put(ImputationJobMapKeys.RUN_SET_ID, runSetStateResponse.getRunSetId());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about submitting a run set
    return StepResult.getStepResultSuccess();
  }
}
