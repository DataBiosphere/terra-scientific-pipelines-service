package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
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

public class PrepareInputsStep implements Step {
  private final PipelinesService pipelinesService;
  private final Logger logger = LoggerFactory.getLogger(PrepareInputsStep.class);

  public PrepareInputsStep(PipelinesService pipelinesService) {
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
        RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS);

    PipelinesEnum pipelineEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    Map<String, Object> userProvidedPipelineInputs =
        inputParameters.get(
            RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS, new TypeReference<>() {});

    Map<String, Object> allPipelineInputs =
        pipelinesService.constructRawInputs(pipelineEnum, userProvidedPipelineInputs);

    List<PipelineInputDefinition> allInputDefinitions =
        pipelinesService.getAllPipelineInputDefinitions(pipelineEnum);

    // storage container stuff
    List<String> keysToPrependWithStorageURL =
        List.of("ref_dict", "reference_panel_path", "genetic_maps_path");
    // in future this will be pulled from the workspace
    String workspaceStorageContainerUri =
        "https://lz8b0d07a4d28c13150a1a12.blob.core.windows.net/sc-94fd136b-4231-4e80-ab0c-76d8a2811066";

    // use input definitions to cast and format all the inputs
    Map<String, Object> formattedPipelineInputs = new HashMap<>();
    for (PipelineInputDefinition inputDefinition : allInputDefinitions) {
      String keyName = inputDefinition.getName();
      PipelineInputTypesEnum pipelineInputType =
          PipelineInputTypesEnum.valueOf(inputDefinition.getType());
      String rawValue;
      if (keysToPrependWithStorageURL.contains(keyName)) {
        rawValue = workspaceStorageContainerUri + allPipelineInputs.get(keyName).toString();
      } else {
        rawValue = allPipelineInputs.get(keyName).toString();
      }
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
