package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.*;
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
  ObjectMapper objectMapper = new ObjectMapper();

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
    Pipeline dbResult = pipelinesRepository.findByName(pipelineName.getValue());
    if (dbResult == null) {
      throw new NotFoundException(
          String.format("Pipeline not found for pipelineName %s", pipelineName));
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

  public void validateInputs(PipelinesEnum pipelineName, Object inputs) {
    Pipeline pipeline = getPipeline(pipelineName);
    List<PipelineInputDefinition> inputDefinitions = pipeline.getPipelineInputDefinitions();

    LinkedHashMap<String, Object> inputsMap = castInputsToMap(inputs);

    ArrayList<String> errorMessages =
        new ArrayList<>(validateRequiredInputs(inputDefinitions, inputsMap));

    errorMessages.addAll(validateInputTypes(inputDefinitions, inputsMap));

    checkForExtraInputs(pipeline, inputDefinitions, inputsMap);

    if (!errorMessages.isEmpty()) {
      throw new ValidationException(
          String.format("Problems with pipelineInputs: %s", String.join("; ", errorMessages)));
    }
  }

  private LinkedHashMap<String, Object> castInputsToMap(Object inputs) {
    try {
      return objectMapper.convertValue(inputs, LinkedHashMap.class);
    } catch (IllegalArgumentException e) {
      throw new ValidationException("pipelineInputs must be a JSON object");
    }
  }

  /**
   * Validate that all required inputs are present in the inputsMap
   *
   * @param inputDefinitions - list of input definitions for a pipeline
   * @param inputsMap - map of inputs to validate
   * @return list of egrror messages for missing required inputs
   */
  public List<String> validateRequiredInputs(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    ArrayList<String> errorMessages = new ArrayList<>();
    inputDefinitions.stream()
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
    ArrayList<String> errorMessages = new ArrayList<>();
    inputDefinitions.forEach(
        inputDefinition -> {
          String inputName = inputDefinition.getName();
          if (inputsMap.containsKey(inputName)) {
            PipelineInputTypesEnum inputType =
                PipelineInputTypesEnum.valueOf(inputDefinition.getType().toUpperCase());
            try {
              inputType.cast(
                  inputName, inputsMap.get(inputName)); // cast method includes a null check
            } catch (ValidationException e) { // custom message from PipelineInputTypesEnum
              errorMessages.add(e.getMessage());
            } catch (IllegalArgumentException e) {
              errorMessages.add(
                  String.format(
                      // note that for security we return the name of the field, not the
                      // user-provided value
                      "%s must be of type %s", inputName, inputDefinition.getType()));
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
      Pipeline pipeline, List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    Set<String> expectedInputNames =
        inputDefinitions.stream().map(PipelineInputDefinition::getName).collect(Collectors.toSet());
    Set<String> providedInputNames = new HashSet<>(inputsMap.keySet());
    providedInputNames.removeAll(expectedInputNames);
    if (!providedInputNames.isEmpty()) {
      logger.warn("Extra inputs provided for pipeline {}: {}", pipeline.getName(), String.join(", ", providedInputNames));
    }
  }
}
