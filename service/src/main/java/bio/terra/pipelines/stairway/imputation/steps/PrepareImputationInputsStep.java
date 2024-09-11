package bio.terra.pipelines.stairway.imputation.steps;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;

import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
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
  private final PipelineRunsService pipelineRunsService;
  private final ImputationConfiguration imputationConfiguration;
  private final Logger logger = LoggerFactory.getLogger(PrepareImputationInputsStep.class);

  public PrepareImputationInputsStep(
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      ImputationConfiguration imputationConfiguration) {
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.imputationConfiguration = imputationConfiguration;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
        RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
        RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME);

    PipelinesEnum pipelineEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    List<PipelineInputDefinition> allInputDefinitions =
        inputParameters.get(
            RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS, new TypeReference<>() {});
    Map<String, Object> userProvidedPipelineInputs =
        inputParameters.get(
            RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS, new TypeReference<>() {});
    String controlWorkspaceStorageContainerUrl =
        inputParameters.get(
            RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME, String.class);
    UUID jobId = UUID.fromString(flightContext.getFlightId());

    Map<String, Object> allPipelineInputs =
        pipelinesService.constructRawInputs(allInputDefinitions, userProvidedPipelineInputs);

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
      String rawValue;
      if (keysToPrependWithStorageURL.contains(keyName)) {
        rawValue = storageWorkspaceStorageContainerUrl + allPipelineInputs.get(keyName).toString();
      } else if (userProvidedInputFileKeys.contains(inputDefinition.getName())) {
        rawValue =
            "%s/%s"
                .formatted(
                    controlWorkspaceStorageContainerUrl,
                    constructDestinationBlobNameForUserInputFile(
                        jobId, allPipelineInputs.get(keyName).toString()));
      } else {
        rawValue = allPipelineInputs.get(keyName).toString();
      }
      // we must cast here, otherwise the inputs will not be properly interpreted later by WDS
      formattedPipelineInputs.put(
          wdlVariableName, pipelineInputType.cast(keyName, rawValue, new TypeReference<>() {}));
    }

    workingMap.put(RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, formattedPipelineInputs);
    logger.info(
        "Constructed and formatted {} pipeline inputs: {}", pipelineEnum, formattedPipelineInputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // this is the first step in RunImputationGcpJobFlight.
    // if undoStep is called it means the flight failed
    // to be moved to a StairwayHook in https://broadworkbench.atlassian.net/browse/TSPS-181

    // set PipelineRun status to FAILED
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID.getKeyName());
    pipelineRunsService.markPipelineRunFailed(
        UUID.fromString(flightContext.getFlightId()),
        inputParameters.get(JobMapKeys.USER_ID.getKeyName(), String.class));

    // increment failed runs counter metric
    PipelinesEnum pipelinesEnum =
        PipelinesEnum.valueOf(
            flightContext
                .getInputParameters()
                .get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    MetricsUtils.incrementPipelineRunFailed(pipelinesEnum);

    return StepResult.getStepResultSuccess();
  }
}
