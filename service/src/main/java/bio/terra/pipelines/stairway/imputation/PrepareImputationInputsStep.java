package bio.terra.pipelines.stairway.imputation;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.RetryException;

/**
 * This step prepares the inputs for the imputation pipeline by assembling the (already validated)
 * user-provided inputs with the service-provided inputs. It adds the storage workspace URL to the
 * service-provided inputs that need it, and then it casts all the inputs according to the type
 * specified in the pipeline input definitions.
 *
 * <p>This step expects the pipeline name and user provided pipeline inputs to be provided in the
 * input parameter map
 *
 * <p>This step constructs the formatted pipeline inputs and stores them in the working map.
 */
public class PrepareImputationInputsStep implements Step {
  private final PipelinesService pipelinesService;
  private final Logger logger = LoggerFactory.getLogger(PrepareImputationInputsStep.class);

  public PrepareImputationInputsStep(PipelinesService pipelinesService) {
    this.pipelinesService = pipelinesService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
        RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS);

    PipelinesEnum pipelineEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    List<PipelineInputDefinition> allInputDefinitions =
        inputParameters.get(
            RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS, new TypeReference<>() {});
    Map<String, Object> userProvidedPipelineInputs =
        inputParameters.get(
            RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS, new TypeReference<>() {});

    Map<String, Object> allPipelineInputs =
        pipelinesService.constructRawInputs(allInputDefinitions, userProvidedPipelineInputs);

    // define input file paths that need to be prepended with the workspace storage URL
    List<String> keysToPrependWithStorageURL =
        List.of("ref_dict", "reference_panel_path", "genetic_maps_path");
    // in future (TSPS-242) this will be generated via WSM from the storage workspace workspace_id
    String workspaceStorageContainerUrl =
        "https://lz8b0d07a4d28c13150a1a12.blob.core.windows.net/sc-94fd136b-4231-4e80-ab0c-76d8a2811066";

    // use input definitions to cast and format all the inputs
    Map<String, Object> formattedPipelineInputs = new HashMap<>();
    for (PipelineInputDefinition inputDefinition : allInputDefinitions) {
      String keyName = inputDefinition.getName();
      PipelineInputTypesEnum pipelineInputType = inputDefinition.getType();
      String rawValue;
      if (keysToPrependWithStorageURL.contains(keyName)) {
        rawValue = workspaceStorageContainerUrl + allPipelineInputs.get(keyName).toString();
      } else {
        rawValue = allPipelineInputs.get(keyName).toString();
      }
      // we must cast here, otherwise the inputs will not be properly interpreted later by WDS
      formattedPipelineInputs.put(
          keyName, pipelineInputType.cast(keyName, rawValue, new TypeReference<>() {}));
    }

    workingMap.put(RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, formattedPipelineInputs);
    logger.info(
        "Constructed and formatted {} pipeline inputs: {}", pipelineEnum, formattedPipelineInputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // no undo for this step
    return StepResult.getStepResultSuccess();
  }
}