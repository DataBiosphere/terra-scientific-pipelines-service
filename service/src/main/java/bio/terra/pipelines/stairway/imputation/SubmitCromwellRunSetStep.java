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

    // launch a cbas submission
    // this is mostly a manually generated run set request definition, we'll want to be able to auto
    // generate this in the future
    RunSetRequest runSetRequest =
        new RunSetRequest()
            .runSetDescription(
                String.format("%s - flight id: %s", pipelineName, flightContext.getFlightId()))
            .runSetName("Flight ID: " + flightContext.getFlightId())
            .methodVersionId(methodVersionId)
            .addWorkflowInputDefinitionsItem(
                new WorkflowInputDefinition()
                    .inputName("HelloWorld.scatter_num")
                    .inputType(
                        new ParameterTypeDefinitionPrimitive()
                            .primitiveType(PrimitiveParameterValueType.INT)
                            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                    .source(
                        new ParameterDefinitionRecordLookup()
                            .recordAttribute("scatter")
                            .type(ParameterDefinition.TypeEnum.RECORD_LOOKUP)))
            .addWorkflowOutputDefinitionsItem(
                new WorkflowOutputDefinition()
                    .outputName("HelloWorld.output_file")
                    .outputType(
                        new ParameterTypeDefinitionArray()
                            .nonEmpty(true)
                            .arrayType(
                                new ParameterTypeDefinitionPrimitive()
                                    .primitiveType(PrimitiveParameterValueType.FILE)
                                    .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                            .type(ParameterTypeDefinition.TypeEnum.ARRAY))
                    .destination(
                        new OutputDestinationRecordUpdate()
                            .recordAttribute("output_file")
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
