package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;
import static bio.terra.pipelines.common.utils.FileUtils.getBlobNameFromTerraWorkspaceStorageUrlGcp;

import bio.terra.common.exception.InternalServerErrorException;
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
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputs;
import bio.terra.rawls.model.Entity;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/* Service to encapsulate the logic for processing pipeline inputs and outputs */
@Service
public class PipelineInputsOutputsService {
  private final GcsService gcsService;

  private final PipelineInputsRepository pipelineInputsRepository;
  private final PipelineOutputsRepository pipelineOutputsRepository;
  private final ObjectMapper objectMapper = new ObjectMapper();

  @Autowired
  public PipelineInputsOutputsService(
      GcsService gcsService,
      PipelineInputsRepository pipelineInputsRepository,
      PipelineOutputsRepository pipelineOutputsRepository) {
    this.gcsService = gcsService;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.pipelineOutputsRepository = pipelineOutputsRepository;
  }

  /**
   * Generate signed PUT urls and curl commands for each user-provided file input in the pipeline.
   *
   * <p>Each user-provided file input (assumed to be a path to a local file) is translated into a
   * write-only (PUT) signed url in a location in the pipeline workspace storage container, in a
   * directory defined by the jobId.
   *
   * <p>This signed url along with the source file path provided by the user are used to generate a
   * curl command that the user can run to upload the file to the location in the pipeline workspace
   * storage container.
   */
  public Map<String, Map<String, String>> prepareFileInputs(
      Pipeline pipeline, UUID jobId, Map<String, Object> userProvidedInputs) {
    // get the list of files that the user needs to upload
    List<String> fileInputNames =
        pipeline.getPipelineInputDefinitions().stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .filter(p -> p.getType().equals(PipelineVariableTypesEnum.FILE))
            .map(PipelineInputDefinition::getName)
            .toList();

    String googleProjectId = pipeline.getWorkspaceGoogleProject();
    String bucketName = pipeline.getWorkspaceStorageContainerName();
    // generate a map where the key is the input name, and the value is a map containing the
    // write-only PUT signed url for the file and the full curl command to upload the file

    Map<String, Map<String, String>> fileInputsMap = new HashMap<>();
    for (String fileInputName : fileInputNames) {
      String fileInputValue = (String) userProvidedInputs.get(fileInputName);
      String objectName = constructDestinationBlobNameForUserInputFile(jobId, fileInputValue);
      String signedUrl =
          gcsService.generatePutObjectSignedUrl(googleProjectId, bucketName, objectName).toString();

      fileInputsMap.put(
          fileInputName,
          Map.of(
              "signedUrl",
              signedUrl,
              "curlCommand",
              "curl -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s'"
                  .formatted(fileInputValue, signedUrl)));
    }

    return fileInputsMap;
  }

  /** Convert pipelineInputs map to string and save to the pipelineInputs table */
  public void savePipelineInputs(Long pipelineRunId, Map<String, Object> pipelineInputs) {
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(pipelineRunId);
    pipelineInput.setInputs(mapToString(pipelineInputs));
    pipelineInputsRepository.save(pipelineInput);
  }

  /** Retrieve pipeline inputs from the pipelineInputs table and convert to a map */
  public Map<String, Object> retrievePipelineInputs(PipelineRun pipelineRun) {
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
  public Map<String, String> extractPipelineOutputsFromEntity(
      List<PipelineOutputDefinition> pipelineOutputDefinitions, Entity entity) {
    Map<String, String> outputs = new HashMap<>();
    for (PipelineOutputDefinition outputDefinition : pipelineOutputDefinitions) {
      String keyName = outputDefinition.getName();
      String wdlVariableName = outputDefinition.getWdlVariableName();
      String outputValue =
          (String)
              entity
                  .getAttributes()
                  .get(wdlVariableName); // .get() returns null if the key is missing, or if the
      // value is empty
      if (outputValue == null) {
        throw new InternalServerErrorException(
            "Output %s is empty or missing".formatted(wdlVariableName));
      }
      outputs.put(keyName, outputValue);
    }
    return outputs;
  }

  /**
   * Extract the pipeline outputs from a pipelineRun object, create signed GET (read-only) urls for
   * each file, and return an ApiPipelineRunOutputs object with the outputs.
   *
   * @param pipelineRun object from the pipelineRunsRepository
   * @return ApiPipelineRunOutputs
   */
  public ApiPipelineRunOutputs formatPipelineRunOutputs(PipelineRun pipelineRun) {
    Map<String, Object> outputsMap =
        stringToMap(
            pipelineOutputsRepository.findPipelineOutputsByJobId(pipelineRun.getId()).getOutputs());

    // currently all outputs are paths that will need a signed url
    String workspaceStorageContainerName = pipelineRun.getWorkspaceStorageContainerName();
    outputsMap.replaceAll(
        (k, v) ->
            gcsService
                .generateGetObjectSignedUrl(
                    pipelineRun.getWorkspaceGoogleProject(),
                    workspaceStorageContainerName,
                    getBlobNameFromTerraWorkspaceStorageUrlGcp(
                        (String) v, workspaceStorageContainerName))
                .toString());

    ApiPipelineRunOutputs apiPipelineRunOutputs = new ApiPipelineRunOutputs();
    apiPipelineRunOutputs.putAll(outputsMap);
    return apiPipelineRunOutputs;
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
