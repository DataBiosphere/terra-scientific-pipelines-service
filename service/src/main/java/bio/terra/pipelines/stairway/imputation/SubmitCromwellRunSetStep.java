package bio.terra.pipelines.stairway.imputation;

import bio.terra.cbas.model.*;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.RetryException;
import java.util.ArrayList;
import java.util.Arrays;
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

  public SubmitCromwellRunSetStep(CbasService cbasService, SamService samService) {
    this.cbasService = cbasService;
    this.samService = samService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    // validate and extract parameters from input map
    FlightMap inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.WDL_METHOD_NAME);

    PipelinesEnum pipelineName =
        inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), PipelinesEnum.class);
    String wdlMethodName = inputParameters.get(RunImputationJobFlightMapKeys.WDL_METHOD_NAME, String.class);

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

    // input definitions - hardcoded for now
    // in future this will be pulled from the workspace
    String workspaceStorageContainerUri = "https://lz8b0d07a4d28c13150a1a12.blob.core.windows.net/sc-94fd136b-4231-4e80-ab0c-76d8a2811066";

    // in future these will be pulled from the working map
    List<String> contigsInputValue =
            new ArrayList<>(
                    Arrays.asList(
                            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                            "chr20", "chr21", "chr22"));
    String refDictInputValue = workspaceStorageContainerUri +
        "/hg38/Homo_sapiens_assembly38.dict";
    String geneticMapsPathInputValue = workspaceStorageContainerUri +
            "/plink-genetic-maps/GRCh38_fixed/";
    String referencePanelPathInputValue = workspaceStorageContainerUri +
            "/hg38/Homo_sapiens_assembly38.dict";

    // launch a cbas submission
    // this is mostly a manually generated run set request definition, we'll want to be able to auto
    // generate this in the future
    RunSetRequest runSetRequest =
        new RunSetRequest()
            .runSetDescription(
                String.format("%s - flight id: %s", pipelineName, flightContext.getFlightId()))
            .runSetName("%s - flightId %s".formatted(wdlMethodName, flightContext.getFlightId()))
            .methodVersionId(methodVersionId)
            // INPUTS
            // contigs input
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("%s.contigs".formatted(wdlMethodName))
                    .inputType(
                        new ParameterTypeDefinitionArray()
                            .nonEmpty(true)
                            .arrayType(
                                new ParameterTypeDefinitionPrimitive()
                                    .primitiveType(PrimitiveParameterValueType.STRING)
                                    .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                            .type(ParameterTypeDefinition.TypeEnum.ARRAY))
                    .source(
                        new ParameterDefinitionLiteralValue()
                            .parameterValue(contigsInputValue)
                            .type(ParameterDefinition.TypeEnum.LITERAL)))
            // genetic_maps_path input
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("%s.genetic_maps_path".formatted(wdlMethodName))
                    .inputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.STRING)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .source(
                        new ParameterDefinitionLiteralValue()
                            .parameterValue(geneticMapsPathInputValue)
                            .type(ParameterDefinition.TypeEnum.LITERAL)))
            // ref_dict input
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("%s.ref_dict".formatted(wdlMethodName))
                    .inputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .source(
                        new ParameterDefinitionLiteralValue()
                            .parameterValue(refDictInputValue)
                            .type(ParameterDefinition.TypeEnum.LITERAL)))
            // reference_panel_path input
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("%s.reference_panel_path".formatted(wdlMethodName))
                    .inputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .source(
                        new ParameterDefinitionLiteralValue()
                            .parameterValue(referencePanelPathInputValue)
                            .type(ParameterDefinition.TypeEnum.LITERAL)))
            // multi_sample_vcf input
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("%s.multi_sample_vcf".formatted(wdlMethodName))
                    .inputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .source(
                        new ParameterDefinitionRecordLookup()
                            .recordAttribute("multi_sample_vcf")
                            .type(ParameterDefinition.TypeEnum.RECORD_LOOKUP)))
            // output_basename input
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("%s.output_basename".formatted(wdlMethodName))
                    .inputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.FILE)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .source(
                        new ParameterDefinitionRecordLookup()
                            .recordAttribute("output_basename")
                            .type(ParameterDefinition.TypeEnum.RECORD_LOOKUP)))
            // OUTPUTS
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
            .wdsRecords(
                new WdsRecordSet()
                    .recordType(pipelineName.getValue())
                    .addRecordIdsItem(flightContext.getFlightId()));
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
