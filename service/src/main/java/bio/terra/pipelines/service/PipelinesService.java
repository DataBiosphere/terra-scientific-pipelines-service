package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.model.PipelineDefinition;
import bio.terra.pipelines.service.pipeline.PipelineDefinitionProvider;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.constraints.NotNull;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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

  private final PipelineDefinitionProvider pipelineDefinitionProvider;
  private final PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;
  private final RawlsService rawlsService;
  private final SamService samService;

  // \\d means digit [0-9]; \\w means word-like [a-zA-Z0-9_]
  private static final String SEM_VER_REGEX_STRING = "^(\\d+|\\w*v\\d+)\\.(\\d+)\\.(\\d+)$";

  @Autowired
  public PipelinesService(
      PipelineDefinitionProvider pipelineDefinitionProvider,
      PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository,
      RawlsService rawlsService,
      SamService samService) {
    this.pipelineDefinitionProvider = pipelineDefinitionProvider;
    this.pipelineRuntimeMetadataRepository = pipelineRuntimeMetadataRepository;
    this.rawlsService = rawlsService;
    this.samService = samService;
  }

  /**
   * Get all pipelines, optionally including hidden pipelines.
   *
   * <p>Pipeline definitions are read from YAML configuration. Runtime metadata (hidden status,
   * workspace details) is merged from the database.
   *
   * @param showHidden - whether the call should return hidden pipelines (this should only happen
   *     for admin users)
   * @return List of pipelines ordered by name (alphabetically) and version (descending)
   */
  public List<Pipeline> getPipelines(boolean showHidden) {
    logger.info("Get all Pipelines");
    List<Pipeline> pipelines = new ArrayList<>();

    for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
      List<PipelineDefinition> definitions =
          pipelineDefinitionProvider.getDefinitionsForPipeline(pipelineEnum);
      for (PipelineDefinition definition : definitions) {
        PipelineRuntimeMetadata runtimeMetadata =
            pipelineRuntimeMetadataRepository.findById(definition.getPipelineKey()).orElse(null);

        if (runtimeMetadata != null && runtimeMetadata.isHidden() && !showHidden) {
          continue;
        }

        pipelines.add(Pipeline.fromDefinitionAndRuntime(definition, runtimeMetadata));
      }
    }

    return pipelines;
  }

  /**
   * Get a specific pipeline by name and version. If version is null, get the latest non hidden
   * version. Will only return non-hidden pipelines unless the user is an admin and a version is
   * specified.
   *
   * <p>Pipeline definitions are read from YAML configuration. Runtime metadata (hidden status,
   * workspace details) is merged from the database.
   *
   * @param pipelineName - name of the pipeline to retrieve
   * @param pipelineVersion - version of the pipeline to retrieve, or null to get the latest version
   * @param showHidden - whether the call should return hidden pipelines (this should only happen
   *     for admin users)
   * @return - the requested Pipeline if it exists
   */
  public Pipeline getPipeline(
      PipelinesEnum pipelineName, Integer pipelineVersion, boolean showHidden) {
    logger.info(
        "Get a specific pipeline for pipelineName {} and version {}",
        pipelineName,
        pipelineVersion);
    if (pipelineVersion == null) {
      return getLatestPipeline(pipelineName, showHidden);
    }

    PipelineDefinition definition =
        pipelineDefinitionProvider.getPipelineDefinition(pipelineName, pipelineVersion);
    PipelineRuntimeMetadata runtimeMetadata =
        pipelineRuntimeMetadataRepository.findById(definition.getPipelineKey()).orElse(null);

    if (runtimeMetadata != null && runtimeMetadata.isHidden() && !showHidden) {
      throw new NotFoundException(
          "Pipeline not found for pipelineName %s and version %s"
              .formatted(pipelineName, pipelineVersion));
    }

    return Pipeline.fromDefinitionAndRuntime(definition, runtimeMetadata);
  }

  /**
   * Get the latest version of a pipeline.
   *
   * @param pipelineName name of the pipeline
   * @param showHidden whether to include hidden pipelines
   * @return the latest pipeline
   */
  private Pipeline getLatestPipeline(PipelinesEnum pipelineName, boolean showHidden) {
    logger.info("Get the latest pipeline for pipelineName {}", pipelineName);
    List<PipelineDefinition> definitions =
        pipelineDefinitionProvider.getDefinitionsForPipeline(pipelineName);
    if (definitions.isEmpty()) {
      throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
    }

    for (PipelineDefinition definition : definitions) {
      PipelineRuntimeMetadata runtimeMetadata =
          pipelineRuntimeMetadataRepository.findById(definition.getPipelineKey()).orElse(null);

      if (runtimeMetadata != null && runtimeMetadata.isHidden() && !showHidden) {
        continue;
      }

      return Pipeline.fromDefinitionAndRuntime(definition, runtimeMetadata);
    }

    throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
  }

  /**
   * Resolve a pipeline by canonical pipelineKey ({pipeline_name}_v{version}).
   *
   * @param pipelineKey canonical pipeline key
   * @return merged pipeline definition/runtime metadata
   */
  public Pipeline getPipelineByKey(String pipelineKey) {
    PipelinesEnum pipelineName = PipelinesEnum.nameFromPipelineKey(pipelineKey);
    int pipelineVersion = PipelinesEnum.versionFromPipelineKey(pipelineKey);
    return getPipeline(pipelineName, pipelineVersion, true);
  }

  /**
   * Helper method to get the latest pipeline without filtering by hidden status. Used for lookup
   * within admin operations.
   *
   * @param pipelineName the pipeline name
   * @return the latest pipeline definition
   */
  public Pipeline getLatestPipeline(PipelinesEnum pipelineName) {
    return getLatestPipeline(pipelineName, true);
  }

  public Pipeline adminUpdatePipelineWorkspace(
      PipelinesEnum pipelineName,
      Integer pipelineVersion,
      Boolean isHidden,
      @NotNull String workspaceBillingProject,
      @NotNull String workspaceName,
      @NotNull String toolVersion) {
    WorkspaceDetails workspaceDetails =
        rawlsService.getWorkspaceDetails(
            samService.getTeaspoonsServiceAccountToken(), workspaceBillingProject, workspaceName);
    String workspaceStorageContainerUrl = rawlsService.getWorkspaceBucketName(workspaceDetails);
    String workspaceGoogleProject = rawlsService.getWorkspaceGoogleProject(workspaceDetails);

    Pipeline pipeline = getPipeline(pipelineName, pipelineVersion, true);

    // ensure toolVersion follows semantic versioning regex
    final Pattern pattern = Pattern.compile(SEM_VER_REGEX_STRING);
    final Matcher matcher = pattern.matcher(toolVersion);
    if (!matcher.matches()) {
      throw new ValidationException(
          String.format(
              "toolVersion %s does not follow semantic versioning regex %s",
              toolVersion, SEM_VER_REGEX_STRING));
    }

    // Build an updated model with the new runtime values
    Pipeline updatedPipeline =
        pipeline.toBuilder()
            .workspaceBillingProject(workspaceBillingProject)
            .workspaceName(workspaceName)
            .workspaceStorageContainerName(workspaceStorageContainerUrl)
            .workspaceGoogleProject(workspaceGoogleProject)
            .toolVersion(toolVersion)
            .hidden(isHidden != null ? isHidden : pipeline.isHidden())
            .build();

    PipelineRuntimeMetadata runtimeMetadata = new PipelineRuntimeMetadata();
    runtimeMetadata.setPipelineKey(updatedPipeline.getPipelineKey());
    runtimeMetadata.setHidden(updatedPipeline.isHidden());
    runtimeMetadata.setToolVersion(updatedPipeline.getToolVersion());
    runtimeMetadata.setWorkspaceBillingProject(updatedPipeline.getWorkspaceBillingProject());
    runtimeMetadata.setWorkspaceName(updatedPipeline.getWorkspaceName());
    runtimeMetadata.setWorkspaceStorageContainerName(
        updatedPipeline.getWorkspaceStorageContainerName());
    runtimeMetadata.setWorkspaceGoogleProject(updatedPipeline.getWorkspaceGoogleProject());
    pipelineRuntimeMetadataRepository.save(runtimeMetadata);

    return updatedPipeline;
  }
}
