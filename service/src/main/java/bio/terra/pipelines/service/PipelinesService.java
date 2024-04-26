package bio.terra.pipelines.service;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import java.util.*;
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
    Pipeline dbResult = pipelinesRepository.findByName(pipelineName.getValue());
    if (dbResult == null) {
      throw new IllegalArgumentException(
          String.format("Pipeline not found for pipelineName %s", pipelineName));
    }
    return dbResult;
  }

  /**
   * This method is meant to only be called by an admin endpoint to update the control workspace id
   * for a pipeline
   *
   * @param pipelineName - nanme of pipeline to update
   * @param workspaceId - UUID of workspace to update to
   * @return
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

    // if no inputs are required, nothing to validate
    if (inputDefinitions.isEmpty()) {
      return;
    }

    LinkedHashMap<String, Object> inputsMap = castInputsToMap(inputs);

    ArrayList<String> errorMessages =
        new ArrayList<>(validateRequiredInputsArePresent(inputDefinitions, inputsMap));

    errorMessages.addAll(validateInputTypes(inputDefinitions, inputsMap));

    if (!errorMessages.isEmpty()) {
      throw new IllegalArgumentException(String.join("; ", errorMessages));
    }
  }

  private LinkedHashMap<String, Object> castInputsToMap(Object inputs) {
    try {
      return (LinkedHashMap<String, Object>) inputs;
    } catch (ClassCastException e) {
      throw new IllegalArgumentException(
          "Pipeline inputs must be in the format {\"input1\": \"value1\", \"input2\": \"value2\"...}");
    }
  }

  /**
   * Validate that all required inputs are present in the inputsMap
   *
   * @param inputDefinitions - list of input definitions for a pipeline
   * @param inputsMap - map of inputs to validate
   * @return list of error messages for missing required inputs
   */
  public List<String> validateRequiredInputsArePresent(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    ArrayList<String> errorMessages = new ArrayList<>();
    inputDefinitions.stream()
        .filter(PipelineInputDefinition::getIsRequired)
        .forEach(
            inputDefinition -> {
              String inputName = inputDefinition.getName();
              if (!(inputsMap.containsKey(inputName)) || inputsMap.get(inputName) == null) {
                errorMessages.add(String.format("pipelineInput %s is required", inputName));
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
          if (inputsMap.containsKey(inputName) && inputsMap.get(inputName) != null) {
            PipelineInputTypesEnum inputType =
                PipelineInputTypesEnum.valueOf(inputDefinition.getType().toUpperCase());
            try {
              inputType.cast(inputsMap.get(inputName));
            } catch (Exception e) {
              errorMessages.add(
                  String.format(
                      // note that for security we return the name of the field, not the
                      // user-provided value
                      "pipelineInput %s must be of type %s", inputName, inputDefinition.getType()));
            }
          }
        });
    return errorMessages;
  }
}
