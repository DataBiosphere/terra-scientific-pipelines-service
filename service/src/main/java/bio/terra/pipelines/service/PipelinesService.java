package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;

import bio.terra.cbas.model.OutputDestination;
import bio.terra.cbas.model.OutputDestinationRecordUpdate;
import bio.terra.cbas.model.ParameterDefinition;
import bio.terra.cbas.model.ParameterDefinitionRecordLookup;
import bio.terra.cbas.model.ParameterTypeDefinition;
import bio.terra.cbas.model.ParameterTypeDefinitionArray;
import bio.terra.cbas.model.ParameterTypeDefinitionPrimitive;
import bio.terra.cbas.model.PrimitiveParameterValueType;
import bio.terra.cbas.model.WorkflowInputDefinition;
import bio.terra.cbas.model.WorkflowOutputDefinition;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.rawls.model.WorkspaceDetails;
import com.fasterxml.jackson.core.type.TypeReference;
import jakarta.validation.constraints.NotNull;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.springframework.validation.annotation.Validated;

/** The Pipelines Service manages information about the service's available Scientific Pipelines. */
@Component
@Validated
public class PipelinesService {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesService.class);

  private final PipelinesRepository pipelinesRepository;
  private final RawlsService rawlsService;
  private final SamService samService;

  private static final String SEM_VER_REGEX_STRING =
      "^(0|[1-9]\\d*|[a-zA-Z_]*v\\d+)\\.(0|[1-9]\\d*)\\.(0|[1-9]\\d*)$";

  @Autowired
  public PipelinesService(
      PipelinesRepository pipelinesRepository, RawlsService rawlsService, SamService samService) {
    this.pipelinesRepository = pipelinesRepository;
    this.rawlsService = rawlsService;
    this.samService = samService;
  }

  public List<Pipeline> getPipelines() {
    logger.info("Get all Pipelines");
    return pipelinesRepository.findAll();
  }

  public Pipeline getPipeline(PipelinesEnum pipelineName, Integer pipelineVersion) {
    logger.info(
        "Get a specific pipeline for pipelineName {} and version {}",
        pipelineName,
        pipelineVersion);
    if (pipelineVersion == null) {
      return getLatestPipeline(pipelineName);
    }
    Pipeline dbResult = pipelinesRepository.findByNameAndVersion(pipelineName, pipelineVersion);
    if (dbResult == null) {
      throw new NotFoundException(
          "Pipeline not found for pipelineName %s and version %s"
              .formatted(pipelineName, pipelineVersion));
    }
    return dbResult;
  }

  public Pipeline getLatestPipeline(PipelinesEnum pipelineName) {
    logger.info("Get the latest pipeline for pipelineName {}", pipelineName);
    Pipeline dbResult = pipelinesRepository.findFirstByNameOrderByVersionDesc(pipelineName);
    if (dbResult == null) {
      throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
    }
    return dbResult;
  }

  public Pipeline getPipelineById(Long pipelineId) {
    logger.info("Get a specific pipeline for pipelineId {}", pipelineId);
    Pipeline dbResult = pipelinesRepository.findById(pipelineId).orElse(null);
    if (dbResult == null) {
      throw new NotFoundException("Pipeline not found for pipelineId %s".formatted(pipelineId));
    }
    return dbResult;
  }

  /**
   * This method is meant to only be called by an admin endpoint to update pipeline parameters such
   * as control workspace information and wdl method version.
   *
   * <p>Calls Rawls to fetch control workspace metadata based on the workspaceBillingProject and
   * workspaceName.
   *
   * @param pipelineName - name of pipeline to update
   * @param workspaceBillingProject - workspace billing project to update to
   * @param workspaceName - workspace name to update to
   * @param wdlMethodVersion - version of wdl expected to run for corresponding pipeline. must align
   *     with pipeline version
   */
  public Pipeline adminUpdatePipelineWorkspace(
      PipelinesEnum pipelineName,
      Integer pipelineVersion,
      @NotNull String workspaceBillingProject,
      @NotNull String workspaceName,
      @NotNull String wdlMethodVersion) {
    WorkspaceDetails workspaceDetails =
        rawlsService.getWorkspaceDetails(
            samService.getTeaspoonsServiceAccountToken(), workspaceBillingProject, workspaceName);
    String workspaceStorageContainerUrl = rawlsService.getWorkspaceBucketName(workspaceDetails);
    String workspaceGoogleProject = rawlsService.getWorkspaceGoogleProject(workspaceDetails);

    Pipeline pipeline = getPipeline(pipelineName, pipelineVersion);
    pipeline.setWorkspaceBillingProject(workspaceBillingProject);
    pipeline.setWorkspaceName(workspaceName);
    pipeline.setWorkspaceStorageContainerName(workspaceStorageContainerUrl);
    pipeline.setWorkspaceGoogleProject(workspaceGoogleProject);

    // ensure wdlMethodVersion follows semantic versioning regex (can be preceded by a string ending
    // in v)
    final Pattern pattern = Pattern.compile(SEM_VER_REGEX_STRING);
    final Matcher matcher = pattern.matcher(wdlMethodVersion);
    if (!matcher.matches()) {
      throw new ValidationException(
          String.format(
              "wdlMethodVersion %s does not follow semantic versioning regex %s",
              wdlMethodVersion, SEM_VER_REGEX_STRING));
    }

    // ensure that major version of wdlMethodVersion matches the value of the pipeline version.
    // split wdlMethodVersion by 'v' and take the last element of the resulting array
    String wdlMethodExtractedSemVer =
        wdlMethodVersion.split("v")[wdlMethodVersion.split("v").length - 1];
    if (pipeline.getVersion().equals(Integer.parseInt(wdlMethodExtractedSemVer.split("\\.")[0]))) {
      pipeline.setWdlMethodVersion(wdlMethodVersion);
    } else {
      throw new ValidationException(
          String.format(
              "wdlMethodVersion %s does not align with pipeline version %s. The major version of wdlMethodVersion must match pipeline version",
              wdlMethodVersion, pipeline.getVersion()));
    }
    pipelinesRepository.save(pipeline);
    return pipeline;
  }

  /**
   * Validate user-provided inputs for a pipeline. Validation includes a check for required inputs
   * and type checks for all inputs. See IMPLEMENTATION_NOTES.md for more details. If any inputs
   * fail validation, a ValidationException is thrown listing all problems. Extra inputs that are
   * not defined in the pipeline are logged at the WARN level.
   *
   * @param allInputDefinitions - all the input definitions for a pipeline
   * @param inputsMap - user-provided inputs Map<String,Object> to validate
   */
  public void validateUserProvidedInputs(
      List<PipelineInputDefinition> allInputDefinitions, Map<String, Object> inputsMap) {
    List<PipelineInputDefinition> userProvidedInputDefinitions =
        extractUserProvidedInputDefinitions(allInputDefinitions);

    List<String> errorMessages =
        new ArrayList<>(validateRequiredInputs(userProvidedInputDefinitions, inputsMap));

    errorMessages.addAll(validateInputTypes(userProvidedInputDefinitions, inputsMap));

    checkForExtraInputs(userProvidedInputDefinitions, inputsMap);

    if (!errorMessages.isEmpty()) {
      throw new ValidationException(
          "Problem(s) with pipelineInputs: %s".formatted(String.join("; ", errorMessages)));
    }
  }

  /**
   * Validate that all required inputs are present in the inputsMap
   *
   * @param inputDefinitions - list of input definitions for a pipeline
   * @param inputsMap - map of inputs to validate
   * @return list of error messages for missing required inputs
   */
  public List<String> validateRequiredInputs(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    ArrayList<String> errorMessages = new ArrayList<>();
    inputDefinitions.stream()
        .filter(PipelineInputDefinition::isUserProvided)
        .filter(PipelineInputDefinition::isRequired)
        .forEach(
            inputDefinition -> {
              String inputName = inputDefinition.getName();
              if (!(inputsMap.containsKey(inputName))) {
                errorMessages.add(String.format("%s is required", inputName));
              }
            });
    return errorMessages;
  }

  /**
   * Validate that all present inputs are the correct type. We do not check for required inputs
   * here.
   *
   * @param inputDefinitions - list of input definitions for a pipeline
   * @param inputsMap - map of inputs to validate
   * @return list of error messages for inputs that are not the correct type
   */
  public List<String> validateInputTypes(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    List<String> errorMessages = new ArrayList<>();
    inputDefinitions.forEach(
        inputDefinition -> {
          String inputName = inputDefinition.getName();
          if (inputsMap.containsKey(inputName)) {
            PipelineVariableTypesEnum inputType = inputDefinition.getType();
            String validationErrorMessage =
                inputType.validate(inputDefinition, inputsMap.get(inputName));
            if (validationErrorMessage != null) {
              errorMessages.add(validationErrorMessage);
            }
          }
        });
    return errorMessages;
  }

  /**
   * Check for extra inputs that are not defined in the pipeline; log a warning if there are
   * unexpected inputs
   */
  public void checkForExtraInputs(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    Set<String> expectedInputNames =
        inputDefinitions.stream().map(PipelineInputDefinition::getName).collect(Collectors.toSet());
    Set<String> providedInputNames = new HashSet<>(inputsMap.keySet());
    providedInputNames.removeAll(expectedInputNames);
    if (!providedInputNames.isEmpty()) {
      String concatenatedInputNames = String.join(", ", providedInputNames);
      logger.warn("Found extra inputs: {}", concatenatedInputNames);
    }
  }

  public List<PipelineInputDefinition> extractUserProvidedInputDefinitions(
      List<PipelineInputDefinition> allInputDefinitions) {
    return allInputDefinitions.stream().filter(PipelineInputDefinition::isUserProvided).toList();
  }

  public List<PipelineInputDefinition> extractServiceProvidedInputDefinitions(
      List<PipelineInputDefinition> allInputDefinitions) {
    return allInputDefinitions.stream()
        .filter(Predicate.not(PipelineInputDefinition::isUserProvided))
        .toList();
  }

  public List<String> extractUserProvidedFileInputNames(
      List<PipelineInputDefinition> inputDefinitions) {
    return inputDefinitions.stream()
        .filter(PipelineInputDefinition::isUserProvided)
        .filter(p -> p.getType().equals(PipelineVariableTypesEnum.FILE))
        .map(PipelineInputDefinition::getName)
        .toList();
  }

  /**
   * Combine the user-provided inputs map with the service-provided inputs to create a map of all
   * the inputs for a pipeline. This does not cast the inputs to the correct type, nor does it
   * format any file inputs with storage container URLs.
   *
   * @param allInputDefinitions - all the input definitions for a pipeline
   * @param userProvidedPipelineInputs - the user-provided inputs
   * @return Map<String, Object> allPipelineInputs - the combined inputs
   */
  public Map<String, Object> constructRawInputs(
      List<PipelineInputDefinition> allInputDefinitions,
      Map<String, Object> userProvidedPipelineInputs) {

    Map<String, Object> allPipelineInputs = new HashMap<>(userProvidedPipelineInputs);

    List<PipelineInputDefinition> serviceProvidedInputDefinitions =
        extractServiceProvidedInputDefinitions(allInputDefinitions);

    // add default values for service-provided inputs to the allPipelineInputs map
    serviceProvidedInputDefinitions.forEach(
        inputDefinition -> {
          String inputName = inputDefinition.getName();
          Object inputValue = inputDefinition.getDefaultValue();
          allPipelineInputs.put(inputName, inputValue); // store the string value; will cast later
        });

    logger.info("All pipeline inputs: {}", allPipelineInputs);

    return allPipelineInputs;
  }

  /**
   * Prepare a list of CBAS WorkflowInputDefinitions using RecordLookup (i.e. reading from WDS) for
   * a given pipeline and WDL method name.
   *
   * @param pipelineInputDefinitions
   * @param wdlMethodName
   * @return
   */
  public List<WorkflowInputDefinition> prepareCbasWorkflowInputRecordLookupDefinitions(
      List<PipelineInputDefinition> pipelineInputDefinitions, String wdlMethodName) {

    return pipelineInputDefinitions.stream()
        .map(
            pipelineInputDefinition -> {
              String inputName = pipelineInputDefinition.getWdlVariableName();
              ParameterTypeDefinition parameterTypeDefinition =
                  mapVariableTypeToCbasParameterType(pipelineInputDefinition.getType());
              return new WorkflowInputDefinition()
                  .inputName("%s.%s".formatted(wdlMethodName, inputName))
                  .inputType(parameterTypeDefinition)
                  .source(
                      new ParameterDefinitionRecordLookup()
                          .recordAttribute(inputName)
                          .type(ParameterDefinition.TypeEnum.RECORD_LOOKUP));
            })
        .toList();
  }

  protected ParameterTypeDefinition mapVariableTypeToCbasParameterType(
      PipelineVariableTypesEnum type) {
    return switch (type) {
      case STRING -> new ParameterTypeDefinitionPrimitive()
          .primitiveType(PrimitiveParameterValueType.STRING)
          .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
      case FILE -> new ParameterTypeDefinitionPrimitive()
          .primitiveType(PrimitiveParameterValueType.FILE)
          .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
      case INTEGER -> new ParameterTypeDefinitionPrimitive()
          .primitiveType(PrimitiveParameterValueType.INT)
          .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
      case STRING_ARRAY -> new ParameterTypeDefinitionArray()
          .nonEmpty(true)
          .arrayType(
              new ParameterTypeDefinitionPrimitive()
                  .primitiveType(PrimitiveParameterValueType.STRING)
                  .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
          .type(ParameterTypeDefinition.TypeEnum.ARRAY);
      case FILE_ARRAY -> new ParameterTypeDefinitionArray()
          .nonEmpty(true)
          .arrayType(
              new ParameterTypeDefinitionPrimitive()
                  .primitiveType(PrimitiveParameterValueType.FILE)
                  .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
          .type(ParameterTypeDefinition.TypeEnum.ARRAY);
    };
  }

  /**
   * Prepare a list of CBAS WorkflowOutputDefinitions for a given pipeline and WDL method name.
   *
   * @param pipelineOutputDefinitions
   * @param wdlMethodName
   * @return
   */
  public List<WorkflowOutputDefinition> prepareCbasWorkflowOutputRecordUpdateDefinitions(
      List<PipelineOutputDefinition> pipelineOutputDefinitions, String wdlMethodName) {
    return pipelineOutputDefinitions.stream()
        .map(
            pipelineOutputDefinition -> {
              String outputName = pipelineOutputDefinition.getWdlVariableName();
              return new WorkflowOutputDefinition()
                  .outputName("%s.%s".formatted(wdlMethodName, outputName))
                  .outputType(
                      mapVariableTypeToCbasParameterType(pipelineOutputDefinition.getType()))
                  .destination(
                      new OutputDestinationRecordUpdate()
                          .recordAttribute(outputName)
                          .type(OutputDestination.TypeEnum.RECORD_UPDATE));
            })
        .toList();
  }

  /**
   * Retrieve and format the pipeline inputs for a pipeline. Combine the user-provided inputs with
   * the service-provided inputs to create a map of all the inputs for the pipeline. Apply the
   * following customizations: - use custom (environment-specific) values for certain
   * service-provided inputs - prepend the storage workspace container URL to the service-provided
   * inputs that need it - prepend the control workspace container URL to the user-provided file
   * inputs
   *
   * @param jobId UUID
   * @param allInputDefinitions List<PipelineInputDefinition>
   * @param userProvidedPipelineInputs Map<String, Object>
   * @param controlWorkspaceContainerUrl String
   * @param inputsWithCustomValues Map<String, Object> from pipeline Configuration
   * @param keysToPrependWithStorageWorkspaceContainerUrl List<String> from pipeline Configuration
   * @param storageWorkspaceContainerUrl String from pipeline Configuration
   * @return formattedPipelineInputs Map<String, Object>
   */
  public Map<String, Object> formatPipelineInputs(
      UUID jobId,
      List<PipelineInputDefinition> allInputDefinitions,
      Map<String, Object> userProvidedPipelineInputs,
      String controlWorkspaceContainerUrl,
      Map<String, Object> inputsWithCustomValues,
      List<String> keysToPrependWithStorageWorkspaceContainerUrl,
      String storageWorkspaceContainerUrl) {

    Map<String, Object> allPipelineInputs =
        constructRawInputs(allInputDefinitions, userProvidedPipelineInputs);
    List<String> userProvidedInputFileKeys = extractUserProvidedFileInputNames(allInputDefinitions);
    Map<String, Object> formattedPipelineInputs = new HashMap<>();

    for (PipelineInputDefinition inputDefinition : allInputDefinitions) {
      String keyName = inputDefinition.getName();
      String wdlVariableName = inputDefinition.getWdlVariableName();
      PipelineVariableTypesEnum pipelineInputType = inputDefinition.getType();

      // this should not happen because we have already validated the inputs
      if (!inputsWithCustomValues.containsKey(keyName)
          && !allPipelineInputs.containsKey(keyName)
          && inputDefinition.isRequired()) {
        logger.error("Required input {} is missing", keyName);
        throw new InternalServerErrorException("Required input %s is missing".formatted(keyName));
      }

      // use custom value if present, otherwise use the value from raw inputs (allPipelineInputs),
      // otherwise use default value
      String rawValue =
          inputsWithCustomValues
              .getOrDefault(
                  keyName,
                  allPipelineInputs.getOrDefault(keyName, inputDefinition.getDefaultValue()))
              .toString();

      String processedValue;

      if (keysToPrependWithStorageWorkspaceContainerUrl.contains(keyName)) {
        processedValue = storageWorkspaceContainerUrl + rawValue;
      } else if (userProvidedInputFileKeys.contains(keyName)) {
        processedValue =
            "%s/%s"
                .formatted(
                    controlWorkspaceContainerUrl,
                    constructDestinationBlobNameForUserInputFile(jobId, rawValue));
      } else {
        processedValue = rawValue;
      }

      // we must cast here, otherwise the inputs will not be properly interpreted later by WDS
      formattedPipelineInputs.put(
          wdlVariableName,
          pipelineInputType.cast(keyName, processedValue, new TypeReference<>() {}));
    }

    return formattedPipelineInputs;
  }
}
