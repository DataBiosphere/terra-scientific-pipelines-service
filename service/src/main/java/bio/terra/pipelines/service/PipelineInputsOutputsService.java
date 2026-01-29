package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;
import static bio.terra.pipelines.common.utils.FileUtils.constructFilePath;
import static bio.terra.pipelines.common.utils.FileUtils.getBlobNameFromTerraWorkspaceStorageUrlGcp;
import static bio.terra.pipelines.common.utils.FileUtils.getFileNameFromFullPath;
import static bio.terra.pipelines.common.utils.FileUtils.isCloudFile;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputSignedUrls;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputs;
import bio.terra.rawls.model.Entity;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
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
import org.springframework.stereotype.Service;

/* Service to encapsulate the logic for processing pipeline inputs and outputs */
@Service
public class PipelineInputsOutputsService {
  private static final Logger logger = LoggerFactory.getLogger(PipelineInputsOutputsService.class);

  private final GcsService gcsService;
  private final PipelineInputsRepository pipelineInputsRepository;
  private final PipelineOutputsRepository pipelineOutputsRepository;
  private final ObjectMapper objectMapper;

  @Autowired
  public PipelineInputsOutputsService(
      GcsService gcsService,
      PipelineInputsRepository pipelineInputsRepository,
      PipelineOutputsRepository pipelineOutputsRepository,
      ObjectMapper objectMapper) {
    this.gcsService = gcsService;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.pipelineOutputsRepository = pipelineOutputsRepository;
    this.objectMapper = objectMapper;
  }

  /**
   * Check whether all user-provided FILE inputs for a pipeline are cloud paths.
   *
   * @param pipeline
   * @param userProvidedInputs
   * @return
   */
  public boolean userProvidedInputsAreCloud(
      Pipeline pipeline, Map<String, Object> userProvidedInputs) {
    List<String> fileInputNames = getFileInputKeys(pipeline);
    for (String fileInputName : fileInputNames) {
      String fileInputValue = (String) userProvidedInputs.get(fileInputName);
      if (!isCloudFile(fileInputValue)) {
        return false;
      }
    }
    return true;
  }

  public List<String> getFileInputKeys(Pipeline pipeline) {
    return pipeline.getPipelineInputDefinitions().stream()
        .filter(PipelineInputDefinition::isUserProvided)
        .filter(p -> p.getType().equals(PipelineVariableTypesEnum.FILE))
        .map(PipelineInputDefinition::getName)
        .toList();
  }

  /**
   * Generate signed PUT/POST urls and curl commands for each user-provided file input in the
   * pipeline, given local file inputs.
   *
   * <p>Each user-provided file input (assumed to be a path to a local file) is translated into a
   * write-only (PUT) signed url in a location in the pipeline workspace storage container, in a
   * directory defined by the jobId. If a resumable upload is requested, a resumable POST signed url
   * is generated for each file input instead.
   *
   * <p>This signed url along with the source file path provided by the user are used to generate a
   * curl command that the user can run to upload the file to the location in the pipeline workspace
   * storage container.
   */
  public Map<String, Map<String, String>> prepareLocalFileInputs(
      Pipeline pipeline,
      UUID jobId,
      Map<String, Object> userProvidedInputs,
      boolean useResumableUploads) {
    List<String> fileInputNames = getFileInputKeys(pipeline);

    String googleProjectId = pipeline.getWorkspaceGoogleProject();
    String bucketName = pipeline.getWorkspaceStorageContainerName();
    // generate a map where the key is the input name, and the value is a map containing the
    // write-only PUT signed url for the file and the full curl command to upload the file

    Map<String, Map<String, String>> fileInputsMap = new HashMap<>();
    for (String fileInputName : fileInputNames) {
      String fileInputValue = (String) userProvidedInputs.get(fileInputName);
      String objectName = constructDestinationBlobNameForUserInputFile(jobId, fileInputValue);
      String signedUrl = getSignedUrl(googleProjectId, bucketName, objectName, useResumableUploads);

      fileInputsMap.put(
          fileInputName,
          Map.of(
              "signedUrl",
              signedUrl,
              "curlCommand",
              getCurlCommand(fileInputValue, signedUrl, useResumableUploads)));
    }

    return fileInputsMap;
  }

  private String getSignedUrl(
      String googleProjectId, String bucketName, String objectName, boolean useResumableUploads) {
    if (useResumableUploads) {
      return gcsService
          .generateResumablePostObjectSignedUrl(googleProjectId, bucketName, objectName)
          .toString();
    } else {
      return gcsService
          .generatePutObjectSignedUrl(googleProjectId, bucketName, objectName)
          .toString();
    }
  }

  private String getCurlCommand(
      String fileInputValue, String signedUrl, boolean useResumableUploads) {
    if (useResumableUploads) {
      return "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s $(curl -s -i -X POST -H 'x-goog-resumable: start' '%s' | grep -i '^Location:' | cut -d' ' -f2- | tr -d '\r') | cat"
          .formatted(fileInputValue, signedUrl);
    } else {
      return "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s' | cat"
          .formatted(fileInputValue, signedUrl);
    }
  }

  /** Convert pipelineInputs map to string and save to the pipelineInputs table */
  public void savePipelineInputs(Long pipelineRunId, Map<String, Object> pipelineInputs) {
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(pipelineRunId);
    pipelineInput.setInputs(mapToString(pipelineInputs));
    pipelineInputsRepository.save(pipelineInput);
  }

  /** Retrieve pipeline inputs from the pipelineInputsRepository and convert to a map */
  public Map<String, Object> retrieveUserProvidedInputs(PipelineRun pipelineRun) {
    PipelineInput pipelineInput =
        pipelineInputsRepository
            .findById(pipelineRun.getId())
            .orElseThrow(
                () ->
                    new InternalServerErrorException(
                        "Pipeline inputs not found for jobId %s"
                            .formatted(pipelineRun.getJobId())));
    return stringToMap(pipelineInput.getInputs());
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

    errorMessages.addAll(checkForExtraInputs(userProvidedInputDefinitions, inputsMap));

    if (!errorMessages.isEmpty()) {
      throw new ValidationException(
          "Problem%s with pipelineInputs: %s"
              .formatted(errorMessages.size() > 1 ? "s" : "", String.join("; ", errorMessages)));
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
   * Check for extra inputs that are not defined in the pipeline; return a list of an error message
   * containing extra inputs found, or an empty list if no extra inputs are found.
   *
   * <p>Even though this method will only ever produce a single error message (if any), we return a
   * list because the method that calls it can handle an empty list better than a single string or
   * null.
   *
   * @param inputDefinitions - list of input definitions for a pipeline
   * @param inputsMap - map of inputs to validate
   * @return list of errorMessage string or empty list if no errors
   */
  public List<String> checkForExtraInputs(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    Set<String> expectedInputNames =
        inputDefinitions.stream().map(PipelineInputDefinition::getName).collect(Collectors.toSet());
    Set<String> providedInputNames = new HashSet<>(inputsMap.keySet());
    providedInputNames.removeAll(expectedInputNames);
    if (!providedInputNames.isEmpty()) {
      return List.of(
          "Found extra input%s (%s)"
              .formatted(
                  providedInputNames.size() > 1 ? "s" : "", String.join(", ", providedInputNames)));
    }
    return List.of();
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

  /**
   * Helper method to add default values for any missing inputs based on the provided list of input
   * definitions.
   *
   * @param inputDefinitions - list of input definitions to use for adding default values
   * @param inputsMap - map of inputs to add default values to
   * @return Map<String, Object> inputsMap - the updated inputs map
   */
  public Map<String, Object> addDefaultValuesForMissingInputs(
      List<PipelineInputDefinition> inputDefinitions, Map<String, Object> inputsMap) {
    Map<String, Object> updatedInputsMap = new HashMap<>(inputsMap);

    inputDefinitions.forEach(
        inputDefinition -> {
          String inputName = inputDefinition.getName();
          if (!updatedInputsMap.containsKey(inputName)) {
            // add default value for missing inputs
            updatedInputsMap.put(inputName, inputDefinition.getDefaultValue());
          }
        });

    return updatedInputsMap;
  }

  /**
   * Add default values for any optional user-provided inputs that were not specified by the user.
   *
   * <p>We expect the user provided inputs to have been validated for all required inputs before
   * calling this method.
   *
   * @param allInputDefinitions - all the input definitions for a pipeline
   * @param userProvidedPipelineInputs - the user-provided inputs
   * @return Map<String, Object> userProvidedPipelineInputs - the updated user-provided inputs
   */
  public Map<String, Object> populateDefaultValuesForMissingOptionalUserInputs(
      List<PipelineInputDefinition> allInputDefinitions,
      Map<String, Object> userProvidedPipelineInputs) {

    return addDefaultValuesForMissingInputs(
        extractUserProvidedInputDefinitions(allInputDefinitions), userProvidedPipelineInputs);
  }

  /**
   * Combine the user-provided inputs map with the service-provided inputs to create a map of all
   * the inputs for a pipeline. Use default values to populate the service-provided inputs. (Note
   * that we expect optional user-provided inputs that weren't specified by the user to have already
   * been populated with default values.)
   *
   * <p>Note that this method does not perform any validation on the inputs. We expect the inputs to
   * have been validated (e.g. all required inputs are present and of the correct type, no extra
   * inputs).
   *
   * <p>This does not cast the inputs to the correct type, nor does it format any file inputs with
   * storage container URLs.
   *
   * @param allInputDefinitions - all the input definitions for a pipeline
   * @param userProvidedPipelineInputs - the user-provided inputs
   * @return Map<String, Object> allPipelineInputs - the combined inputs
   */
  Map<String, Object> addServiceProvidedInputs(
      List<PipelineInputDefinition> allInputDefinitions,
      Map<String, Object> userProvidedPipelineInputs) {

    return addDefaultValuesForMissingInputs(
        extractServiceProvidedInputDefinitions(allInputDefinitions), userProvidedPipelineInputs);
  }

  /**
   * Format the pipeline inputs for a pipeline. Apply the following manipulations:
   *
   * <ul>
   *   <li>use custom (environment-specific) values for certain service-provided inputs
   *   <li>prepend the storage workspace container URL to the service-provided inputs that need it
   *   <li>prepend the control workspace container URL to the user-provided file inputs
   *   <li>cast all the inputs according to the type specified in the pipeline input definitions
   * </ul>
   */
  Map<String, Object> formatPipelineInputs(
      Map<String, Object> allRawInputs,
      List<PipelineInputDefinition> allInputDefinitions,
      UUID jobId,
      String controlWorkspaceContainerUrl,
      Map<String, String> inputsWithCustomValues,
      List<String> keysToPrependWithStorageWorkspaceContainerUrl,
      String storageWorkspaceContainerUrl) {
    Map<String, Object> formattedPipelineInputs = new HashMap<>();

    for (PipelineInputDefinition inputDefinition : allInputDefinitions) {
      String keyName = inputDefinition.getName();
      String wdlVariableName = inputDefinition.getWdlVariableName();
      PipelineVariableTypesEnum pipelineInputType = inputDefinition.getType();

      // use custom value if present, otherwise use the value from raw inputs (allRawInputs)
      String rawOrCustomValue =
          (inputsWithCustomValues.containsKey(keyName))
              ? inputsWithCustomValues.get(keyName)
              : allRawInputs.get(keyName).toString();
      String processedValue;

      if (keysToPrependWithStorageWorkspaceContainerUrl.contains(keyName)) {
        // the rawOrCustomValue for this field should start with a / so we don't need to add one
        // here
        processedValue = constructFilePath(storageWorkspaceContainerUrl, rawOrCustomValue);
      } else if (inputDefinition.isUserProvided()
          && inputDefinition.getType().equals(PipelineVariableTypesEnum.FILE)
          && !isCloudFile(rawOrCustomValue)) {
        // user-provided file inputs are formatted with control workspace container url and a custom
        // path
        processedValue =
            constructFilePath(
                controlWorkspaceContainerUrl,
                constructDestinationBlobNameForUserInputFile(jobId, rawOrCustomValue));
      } else {
        processedValue = rawOrCustomValue;
      }

      // we must cast here, otherwise the inputs will not be properly interpreted later by WDS
      formattedPipelineInputs.put(
          wdlVariableName,
          pipelineInputType.cast(keyName, processedValue, new TypeReference<>() {}));
    }

    logger.info("Formatted pipeline inputs: {}", formattedPipelineInputs);

    return formattedPipelineInputs;
  }

  /**
   * Gather and format the inputs for a pipeline. Gathers the user-provided inputs and combines with
   * service-provided inputs, then formats all inputs with any custom values and bucket paths, and
   * casts final values to the appropriate type as defined by input definitions. Uses the
   * wdlVariableName as the key for the formatted inputs so that these values can be used as
   * pipeline inputs.
   *
   * @param jobId UUID
   * @param allInputDefinitions List<PipelineInputDefinition>
   * @param userProvidedPipelineInputs Map<String, Object>
   * @param controlWorkspaceContainerUrl String
   * @param inputsWithCustomValues Map<String, String> from pipeline Configuration
   * @param keysToPrependWithStorageWorkspaceContainerUrl List<String> from pipeline Configuration
   * @param storageWorkspaceContainerUrl String from pipeline Configuration
   * @return formattedPipelineInputs Map<String, Object>
   */
  public Map<String, Object> gatherAndFormatPipelineInputs(
      UUID jobId,
      List<PipelineInputDefinition> allInputDefinitions,
      Map<String, Object> userProvidedPipelineInputs,
      String controlWorkspaceContainerUrl,
      Map<String, String> inputsWithCustomValues,
      List<String> keysToPrependWithStorageWorkspaceContainerUrl,
      String storageWorkspaceContainerUrl) {

    Map<String, Object> allRawInputs =
        addServiceProvidedInputs(allInputDefinitions, userProvidedPipelineInputs);
    return formatPipelineInputs(
        allRawInputs,
        allInputDefinitions,
        jobId,
        controlWorkspaceContainerUrl,
        inputsWithCustomValues,
        keysToPrependWithStorageWorkspaceContainerUrl,
        storageWorkspaceContainerUrl);
  }

  // methods to interact with and format pipeline run outputs

  /**
   * Extract pipeline outputs from a Rawls entity object, converting wdlVariableName (typically
   * snake_case) to outputName (camelCase). Throw an error if any outputs are missing from the
   * entity or empty.
   *
   * @param pipelineOutputDefinitions
   * @param entity
   * @return a map of pipeline outputs
   */
  public Map<String, Object> extractPipelineOutputsFromEntity(
      List<PipelineOutputDefinition> pipelineOutputDefinitions, Entity entity) {
    Map<String, Object> outputs = new HashMap<>();
    for (PipelineOutputDefinition outputDefinition : pipelineOutputDefinitions) {
      String keyName = outputDefinition.getName();
      String wdlVariableName = outputDefinition.getWdlVariableName();
      PipelineVariableTypesEnum outputType = outputDefinition.getType();
      boolean isRequired = outputDefinition.isRequired();
      Object outputValue =
          entity
              .getAttributes()
              .get(wdlVariableName); // .get() returns null if the key is missing, or if the
      // value is empty
      if (isRequired && outputValue == null) {
        throw new InternalServerErrorException(
            "Output %s is empty or missing".formatted(wdlVariableName));
      }
      outputs.put(keyName, outputType.cast(keyName, outputValue, new TypeReference<>() {}));
    }
    return outputs;
  }

  /**
   * Retrieve the pipeline outputs from a pipelineRun object and return an ApiPipelineRunOutputs
   * object containing the outputs with files reduced to file names.
   *
   * <p>We expect the pipeline run to have been confirmed as SUCCEEDED before this is called.
   *
   * @param pipelineRun object from the pipelineRunsRepository
   * @return ApiPipelineRunOutputs
   */
  public ApiPipelineRunOutputs getPipelineRunOutputs(PipelineRun pipelineRun) {
    Map<String, Object> outputsMap =
        stringToMap(
            pipelineOutputsRepository.findPipelineOutputsByJobId(pipelineRun.getId()).getOutputs());

    // for any outputs that are file paths, reduce to just the file name
    Set<String> fileOutputNames = getFileOutputKeys(pipelineRun.getPipeline());

    outputsMap.replaceAll(
        (key, value) ->
            fileOutputNames.contains(key) ? getFileNameFromFullPath((String) value) : value);

    ApiPipelineRunOutputs apiPipelineRunOutputs = new ApiPipelineRunOutputs();
    apiPipelineRunOutputs.putAll(outputsMap);
    return apiPipelineRunOutputs;
  }

  /**
   * Extract the pipeline FILE outputs from a pipelineRun object, create signed GET (read-only) urls
   * for each output file, and return an ApiPipelineRunOutputSignedUrls object containing the signed
   * urls.
   *
   * @param pipelineRun object from the pipelineRunsRepository
   * @return ApiPipelineRunOutputSignedUrls containing signed urls for file outputs
   */
  public ApiPipelineRunOutputSignedUrls generatePipelineRunOutputSignedUrls(
      PipelineRun pipelineRun) {
    Map<String, Object> outputsMap =
        stringToMap(
            pipelineOutputsRepository.findPipelineOutputsByJobId(pipelineRun.getId()).getOutputs());
    Map<String, String> signedUrls = new HashMap<>();

    String workspaceStorageContainerName = pipelineRun.getWorkspaceStorageContainerName();
    // populate signedUrls with signed URLs for each file output
    for (String outputName : getFileOutputKeys(pipelineRun.getPipeline())) {
      String filePath = (String) outputsMap.get(outputName);
      String signedUrl =
          gcsService
              .generateGetObjectSignedUrl(
                  pipelineRun.getWorkspaceGoogleProject(),
                  workspaceStorageContainerName,
                  getBlobNameFromTerraWorkspaceStorageUrlGcp(
                      filePath, workspaceStorageContainerName))
              .toString();
      signedUrls.put(outputName, signedUrl);
    }

    ApiPipelineRunOutputSignedUrls apiPipelineRunOutputWithSignedUrls =
        new ApiPipelineRunOutputSignedUrls();
    apiPipelineRunOutputWithSignedUrls.putAll(signedUrls);
    return apiPipelineRunOutputWithSignedUrls;
  }

  /**
   * Get the set of a pipeline's output definition keys that are of type FILE
   *
   * @param pipeline
   * @return Set<String> of FILE-type output keys
   */
  private Set<String> getFileOutputKeys(Pipeline pipeline) {
    return pipeline.getPipelineOutputDefinitions().stream()
        .filter(def -> def.getType().equals(PipelineVariableTypesEnum.FILE))
        .map(PipelineOutputDefinition::getName)
        .collect(Collectors.toSet());
  }

  /** Convert pipelineOutputs map to string and save to the pipelineOutputs table */
  public void savePipelineOutputs(Long pipelineRunId, Map<String, String> pipelineOutputs) {
    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(pipelineRunId);
    pipelineOutput.setOutputs(mapToString(pipelineOutputs));
    pipelineOutputsRepository.save(pipelineOutput);
  }

  public String mapToString(Map<String, ?> outputsMap) {
    try {
      return objectMapper.writeValueAsString(outputsMap);
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error converting map to string", e);
    }
  }

  public Map<String, Object> stringToMap(String outputsString) {
    try {
      return objectMapper.readValue(outputsString, new TypeReference<>() {});
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error converting string to map", e);
    }
  }
}
