package bio.terra.pipelines.stairway.imputation;

import bio.terra.cbas.model.*;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import java.util.List;
import java.util.UUID;

/**
 * This step submits a run set to cromwell. It first finds the correct method version id to run
 * based on the method name passed to it. It then generates a run set request (currently hardcoded
 * but should be generated in the future) that links that request to wds record(s). It then submits
 * that request to cbas
 *
 * <p>this step expects the cbas uri to be passed in through the working map
 *
 * <p>this step writes run set id to the working map
 */
public class SubmitCromwellRunSetStep implements Step {
  private final CbasService cbasService;
  private final SamService samService;
  private final PipelinesService pipelinesService;

  public SubmitCromwellRunSetStep(
      CbasService cbasService, SamService samService, PipelinesService pipelinesService) {
    this.cbasService = cbasService;
    this.samService = samService;
    this.pipelinesService = pipelinesService;
  }

  @Override
  @SuppressWarnings("java:S2259") // suppress warning for possible NPE - we do validate not null in
  // `validateRequiredEntries`
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.DESCRIPTION.getKeyName(),
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.WDL_METHOD_NAME);

    String description = inputParameters.get(JobMapKeys.DESCRIPTION.getKeyName(), String.class);
    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);
    String wdlMethodName =
        inputParameters.get(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, String.class);

    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, RunImputationJobFlightMapKeys.CBAS_URI);

    String cbasUri = workingMap.get(RunImputationJobFlightMapKeys.CBAS_URI, String.class);

    // grab methodVersionId needed to submit a submission
    MethodListResponse methodListResponse =
        cbasService.getAllMethods(cbasUri, samService.getTspsServiceAccountToken());
    UUID methodVersionId =
        CbasService.getMethodVersionIdFromMethodListResponse(methodListResponse, wdlMethodName);
    if (methodVersionId == null) {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException(
              "couldn't find method version id for method name " + wdlMethodName));
    }

    // prepare cbas submission (run set) object
    // this is mostly a manually generated run set request definition, we'll want to be able to auto
    // generate this in the future
    RunSetRequest runSetRequest =
        new RunSetRequest()
            .runSetDescription(
                "%s (%s) - flight id: %s; description: %s"
                    .formatted(
                        pipelineName, wdlMethodName, flightContext.getFlightId(), description))
            .runSetName("%s - flightId %s".formatted(wdlMethodName, flightContext.getFlightId()))
            .methodVersionId(methodVersionId)
            // OUTPUTS
            // imputed_multi_sample_vcf output
            .addWorkflowOutputDefinitionsItem(
                new WorkflowOutputDefinition()
                    .outputName("%s.imputed_multi_sample_vcf".formatted(wdlMethodName))
                    .outputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .destination(
                        new OutputDestinationRecordUpdate()
                            .recordAttribute("imputed_multi_sample_vcf")
                            .type(OutputDestination.TypeEnum.RECORD_UPDATE)))
            // imputed_multi_sample_vcf_index output
            .addWorkflowOutputDefinitionsItem(
                new WorkflowOutputDefinition()
                    .outputName("%s.imputed_multi_sample_vcf_index".formatted(wdlMethodName))
                    .outputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .destination(
                        new OutputDestinationRecordUpdate()
                            .recordAttribute("imputed_multi_sample_vcf_index")
                            .type(OutputDestination.TypeEnum.RECORD_UPDATE)))
            // chunks_info output
            .addWorkflowOutputDefinitionsItem(
                new WorkflowOutputDefinition()
                    .outputName("%s.chunks_info".formatted(wdlMethodName))
                    .outputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .destination(
                        new OutputDestinationRecordUpdate()
                            .recordAttribute("chunks_info")
                            .type(OutputDestination.TypeEnum.RECORD_UPDATE)))
            // failed_chunks output
            .addWorkflowOutputDefinitionsItem(
                new WorkflowOutputDefinition()
                    .outputName("%s.failed_chunks".formatted(wdlMethodName))
                    .outputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .destination(
                        new OutputDestinationRecordUpdate()
                            .recordAttribute("failed_chunks")
                            .type(OutputDestination.TypeEnum.RECORD_UPDATE)))
            // n_failed_chunks output
            .addWorkflowOutputDefinitionsItem(
                new WorkflowOutputDefinition()
                    .outputName("%s.n_failed_chunks".formatted(wdlMethodName))
                    .outputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.INT)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .destination(
                        new OutputDestinationRecordUpdate()
                            .recordAttribute("n_failed_chunks")
                            .type(OutputDestination.TypeEnum.RECORD_UPDATE)))
            .wdsRecords(new WdsRecordSet().recordType(null))
            .wdsRecords(
                new WdsRecordSet()
                    .recordType(pipelineName.getValue())
                    .addRecordIdsItem(flightContext.getFlightId()));

    // add inputs
    List<WorkflowInputDefinition> cbasWorkflowInputDefinitions =
        pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(
            pipelinesService.getAllPipelineInputDefinitions(pipelineName), wdlMethodName);
    for (WorkflowInputDefinition workflowInputDefinition : cbasWorkflowInputDefinitions) {
      runSetRequest.addWorkflowInputDefinitionsItem(workflowInputDefinition);
    }

    // launch the submission
    RunSetStateResponse runSetStateResponse =
        cbasService.createRunSet(cbasUri, samService.getTspsServiceAccountToken(), runSetRequest);
    workingMap.put(RunImputationJobFlightMapKeys.RUN_SET_ID, runSetStateResponse.getRunSetId());
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext context) throws InterruptedException {
    // nothing to undo; there's nothing to undo about submitting a run set
    return StepResult.getStepResultSuccess();
  }
}
