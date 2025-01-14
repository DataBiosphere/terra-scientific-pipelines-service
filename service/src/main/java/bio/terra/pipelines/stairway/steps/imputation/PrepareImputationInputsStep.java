package bio.terra.pipelines.stairway.steps.imputation;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
  private final ImputationConfiguration imputationConfiguration;
  private final Logger logger = LoggerFactory.getLogger(PrepareImputationInputsStep.class);

  public PrepareImputationInputsStep(
      PipelinesService pipelinesService, ImputationConfiguration imputationConfiguration) {
    this.pipelinesService = pipelinesService;
    this.imputationConfiguration = imputationConfiguration;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME,
        ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS,
        ImputationJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL);

    PipelinesEnum pipelineEnum =
        PipelinesEnum.valueOf(inputParameters.get(JobMapKeys.PIPELINE_NAME, String.class));
    List<PipelineInputDefinition> allInputDefinitions =
        inputParameters.get(
            ImputationJobMapKeys.PIPELINE_INPUT_DEFINITIONS, new TypeReference<>() {});
    Map<String, Object> userProvidedPipelineInputs =
        inputParameters.get(
            ImputationJobMapKeys.USER_PROVIDED_PIPELINE_INPUTS, new TypeReference<>() {});
    String controlWorkspaceStorageContainerName =
        inputParameters.get(
            ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME, String.class);
    String controlWorkspaceStorageContainerProtocol =
        inputParameters.get(
            ImputationJobMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL, String.class);
    UUID jobId = UUID.fromString(flightContext.getFlightId());

    // construct the control workspace storage URL
    String controlWorkspaceStorageContainerUrl =
        "%s%s"
            .formatted(
                controlWorkspaceStorageContainerProtocol, controlWorkspaceStorageContainerName);

    Map<String, Object> allPipelineInputs =
        pipelinesService.constructRawInputs(allInputDefinitions, userProvidedPipelineInputs);

    // define input keys that have custom values to be read from the config
    Map<String, Object> inputsWithCustomValues =
        imputationConfiguration.getInputsWithCustomValues();

    // define input file paths that need to be prepended with the storage workspace storage URL
    List<String> keysToPrependWithStorageURL =
        imputationConfiguration.getInputKeysToPrependWithStorageUrl();
    String storageWorkspaceStorageContainerUrl =
        imputationConfiguration.getStorageWorkspaceStorageUrl();

    // define the user-provided inputs that need to be prepended with the control workspace storage
    // URL
    List<String> userProvidedInputFileKeys =
        pipelinesService.extractUserProvidedFileInputNames(allInputDefinitions);

    // use input definitions to cast and format all the inputs
    Map<String, Object> formattedPipelineInputs = new HashMap<>();
    for (PipelineInputDefinition inputDefinition : allInputDefinitions) {
      String keyName = inputDefinition.getName();
      String wdlVariableName = inputDefinition.getWdlVariableName();
      PipelineVariableTypesEnum pipelineInputType = inputDefinition.getType();
      String rawValue = allPipelineInputs.get(keyName).toString();
      String processedValue;

      // overwrite rawValue with custom values from config
      if (inputsWithCustomValues.containsKey(keyName)) {
        rawValue = inputsWithCustomValues.get(keyName).toString();
      }

      if (keysToPrependWithStorageURL.contains(keyName)) {
        processedValue = storageWorkspaceStorageContainerUrl + rawValue;
      } else if (userProvidedInputFileKeys.contains(inputDefinition.getName())) {
        processedValue =
            "%s/%s"
                .formatted(
                    controlWorkspaceStorageContainerUrl,
                    constructDestinationBlobNameForUserInputFile(jobId, rawValue));
      } else {
        processedValue = rawValue;
      }
      // we must cast here, otherwise the inputs will not be properly interpreted later by WDS
      formattedPipelineInputs.put(
          wdlVariableName,
          pipelineInputType.cast(keyName, processedValue, new TypeReference<>() {}));
    }

    workingMap.put(ImputationJobMapKeys.ALL_PIPELINE_INPUTS, formattedPipelineInputs);
    logger.info(
        "Constructed and formatted {} pipeline inputs: {}", pipelineEnum, formattedPipelineInputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    return StepResult.getStepResultSuccess();
  }
}
