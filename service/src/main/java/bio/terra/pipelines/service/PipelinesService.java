package bio.terra.pipelines.service;

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
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/** The Pipelines Service manages information about the service's available Scientific Pipelines. */
@Component
public class PipelinesService {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesService.class);

  private final PipelinesRepository pipelinesRepository;

  @Autowired
  public PipelinesService(PipelinesRepository pipelinesRepository) {
    this.pipelinesRepository = pipelinesRepository;
  }

  public List<Pipeline> getPipelines() {
    logger.info("Get all Pipelines");
    return pipelinesRepository.findAll();
  }

  public Pipeline getPipeline(PipelinesEnum pipelineName) {
    logger.info("Get a specific pipeline for pipelineName {}", pipelineName);
    Pipeline dbResult = pipelinesRepository.findByName(pipelineName);
    if (dbResult == null) {
      throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
    }
    return dbResult;
  }

  /**
   * This method is meant to only be called by an admin endpoint to update the control workspace id
   * for a pipeline
   *
   * @param pipelineName - name of pipeline to update
   * @param workspaceId - UUID of workspace to update to
   * @return pipeline with updated workspaceId
   */
  public Pipeline updatePipelineWorkspaceId(PipelinesEnum pipelineName, UUID workspaceId) {
    Pipeline pipeline = getPipeline(pipelineName);
    pipeline.setWorkspaceId(workspaceId);
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
        .filter(PipelineInputDefinition::getUserProvided)
        .filter(PipelineInputDefinition::getIsRequired)
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
            PipelineInputTypesEnum inputType = inputDefinition.getType();
            String validationErrorMessage = inputType.validate(inputName, inputsMap.get(inputName));
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
    return allInputDefinitions.stream().filter(PipelineInputDefinition::getUserProvided).toList();
  }

  public List<PipelineInputDefinition> extractServiceProvidedInputDefinitions(
      List<PipelineInputDefinition> allInputDefinitions) {
    return allInputDefinitions.stream()
        .filter(Predicate.not(PipelineInputDefinition::getUserProvided))
        .toList();
  }

  public List<String> extractUserProvidedFileInputNames(
      List<PipelineInputDefinition> inputDefinitions) {
    return inputDefinitions.stream()
        .filter(PipelineInputDefinition::getUserProvided)
        // note VCF is our only file type currently, and this method doesn't yet support
        // VCF_ARRAY
        .filter(p -> p.getType().equals(PipelineInputTypesEnum.VCF))
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
                  mapInputTypeToCbasParameterType(pipelineInputDefinition.getType());
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

  protected ParameterTypeDefinition mapInputTypeToCbasParameterType(PipelineInputTypesEnum type) {
    return switch (type) {
      case STRING, VCF -> new ParameterTypeDefinitionPrimitive()
          .primitiveType(PrimitiveParameterValueType.STRING)
          .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
      case INTEGER -> new ParameterTypeDefinitionPrimitive()
          .primitiveType(PrimitiveParameterValueType.INT)
          .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
      case STRING_ARRAY, VCF_ARRAY -> new ParameterTypeDefinitionArray()
          .nonEmpty(true)
          .arrayType(
              new ParameterTypeDefinitionPrimitive()
                  .primitiveType(PrimitiveParameterValueType.STRING)
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
                      new ParameterTypeDefinitionPrimitive()
                          .primitiveType(PrimitiveParameterValueType.FILE)
                          .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                  .destination(
                      new OutputDestinationRecordUpdate()
                          .recordAttribute(outputName)
                          .type(OutputDestination.TypeEnum.RECORD_UPDATE));
            })
        .toList();
  }
}
