package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.enumFromPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.versionFromPipelineKey;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.PipelineInputDefinitionConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.PipelineOutputDefinitionConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.WdlBasedPipelineConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.constraints.NotNull;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.springframework.validation.annotation.Validated;

/** The Pipelines Service manages information about the service's Scientific Pipelines. */
@Component
@Validated
public class PipelinesService {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesService.class);

  private final PipelineConfigurations pipelineConfigurations;
  private final PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;
  private final RawlsService rawlsService;
  private final SamService samService;

  // \\d means digit [0-9]; \\w means word-like [a-zA-Z0-9_]
  private static final String SEM_VER_REGEX_STRING = "^(\\d+|\\w*v\\d+)\\.(\\d+)\\.(\\d+)$";

  @Autowired
  public PipelinesService(
      PipelineConfigurations pipelineConfigurations,
      PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository,
      RawlsService rawlsService,
      SamService samService) {
    this.pipelineConfigurations = pipelineConfigurations;
    this.pipelineRuntimeMetadataRepository = pipelineRuntimeMetadataRepository;
    this.rawlsService = rawlsService;
    this.samService = samService;
  }

  /**
   * Get all pipelines, optionally including hidden pipelines.
   *
   * <p>Runtime metadata (hidden status, workspace details) is read from the database to select
   * which pipelines (via pipeline_key) to return. The relevant Pipeline definitions are read from
   * YAML configuration.
   *
   * <p>Pipelines are returned sorted by name (alphabetically) and version (descending).
   *
   * <p>Note that uninitialized pipelines - those without a row in the pipelineRuntimeMetadata table
   * - will not show up in this list, even when called by an admin.
   *
   * @param showHidden - whether the call should return hidden pipelines (this should only happen
   *     for admin users)
   * @return List of pipelines ordered by name (alphabetically) and version (descending)
   */
  public List<Pipeline> getPipelines(boolean showHidden) {
    logger.info("Get all Pipelines");

    List<PipelineRuntimeMetadata> metadataList =
        showHidden
            ? pipelineRuntimeMetadataRepository.findAll()
            : pipelineRuntimeMetadataRepository.findAllByHidden(false);
    // TODO either revert this to getting yaml info first, or make separate adminGetPipelines method
    return metadataList.stream()
        .map(
            metadata -> {
              PipelinesEnum pipelineName = enumFromPipelineKey(metadata.getPipelineKey());
              int version = versionFromPipelineKey(metadata.getPipelineKey());
              WdlBasedPipelineConfiguration config =
                  pipelineConfigurations.getPipelineConfiguration(metadata.getPipelineKey());
              return pipelineFromConfigAndMetadata(pipelineName, version, config, metadata);
            })
        .sorted(
            Comparator.comparing((Pipeline p) -> p.getName().getLowerCaseValue())
                .thenComparing(Comparator.comparingInt(Pipeline::getVersion).reversed()))
        .toList();
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
      return getLatestPipeline(pipelineName);
    }

    PipelineRuntimeMetadata runtimeMetadata =
        pipelineRuntimeMetadataRepository
            .findById(buildPipelineKey(pipelineName, pipelineVersion))
            .orElse(null);

    // admins (i.e. showHidden = true) can get pipelines even if no runtimeMetadata is specified
    if (!showHidden && (runtimeMetadata == null || runtimeMetadata.isHidden())) {
      throw new NotFoundException(
          "Pipeline not found for pipelineName %s and version %s"
              .formatted(pipelineName, pipelineVersion));
    }

    WdlBasedPipelineConfiguration config =
        pipelineConfigurations.getPipelineConfiguration(
            buildPipelineKey(pipelineName, pipelineVersion));

    return pipelineFromConfigAndMetadata(pipelineName, pipelineVersion, config, runtimeMetadata);
  }

  /**
   * Get the latest non-hidden version of a pipeline.
   *
   * @param pipelineName name of the pipeline
   * @return the latest pipeline
   */
  public Pipeline getLatestPipeline(PipelinesEnum pipelineName) {
    logger.info("Get the latest visible pipeline for pipelineName {}", pipelineName);

    // Query runtime metadata for all versions of this pipeline
    String pipelineKeyPrefix = pipelineName.getLowerCaseValue();
    PipelineRuntimeMetadata latestMetadata =
        pipelineRuntimeMetadataRepository
            .findFirstByPipelineKeyStartingWithAndHiddenIsFalseOrderByPipelineKeyDesc(
                pipelineKeyPrefix);

    if (latestMetadata == null) {
      throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
    }

    // Extract version from the pipeline key and fetch the definition
    int latestVersion = versionFromPipelineKey(latestMetadata.getPipelineKey());
    WdlBasedPipelineConfiguration config =
        pipelineConfigurations.getPipelineConfiguration(latestMetadata.getPipelineKey());

    return pipelineFromConfigAndMetadata(pipelineName, latestVersion, config, latestMetadata);
  }

  /**
   * Helper method to construct a Pipeline domain model by merging YAML configuration with database
   * runtime metadata.
   */
  private Pipeline pipelineFromConfigAndMetadata(
      PipelinesEnum name,
      Integer version,
      WdlBasedPipelineConfiguration config,
      PipelineRuntimeMetadata runtimeMetadata) {

    Pipeline.PipelineBuilder builder =
        Pipeline.builder()
            .name(name)
            .version(version)
            .pipelineKey(buildPipelineKey(name, version))
            .displayName(config.getDisplayName())
            .description(config.getDescription())
            .pipelineType(config.getPipelineType())
            .toolName(config.getToolName())
            .inputDefinitions(inputDefinitionsFromConfig(config.getInputDefinitions()))
            .outputDefinitions(outputDefinitionsFromConfig(config.getOutputDefinitions()))
            .memoryRetryMultiplier(config.getMemoryRetryMultiplier())
            .hidden(true);

    if (runtimeMetadata != null) {
      builder
          .hidden(runtimeMetadata.isHidden())
          .toolVersion(runtimeMetadata.getToolVersion())
          .workspaceBillingProject(runtimeMetadata.getWorkspaceBillingProject())
          .workspaceName(runtimeMetadata.getWorkspaceName())
          .workspaceStorageContainerName(runtimeMetadata.getWorkspaceStorageContainerName())
          .workspaceGoogleProject(runtimeMetadata.getWorkspaceGoogleProject())
          .updated(runtimeMetadata.getUpdated());
    }

    return builder.build();
  }

  private List<PipelineInputDefinition> inputDefinitionsFromConfig(
      List<PipelineInputDefinitionConfiguration> inputConfigs) {
    if (inputConfigs == null) {
      return java.util.Collections.emptyList();
    }
    return inputConfigs.stream()
        .map(
            config ->
                PipelineInputDefinition.builder()
                    .name(config.getName())
                    .wdlVariableName(config.getWdlVariableName())
                    .displayName(config.getDisplayName())
                    .description(config.getDescription())
                    .type(config.getType())
                    .isRequired(config.getIsRequired() != null && config.getIsRequired())
                    .userProvided(config.getUserProvided() != null && config.getUserProvided())
                    .defaultValue(config.getDefaultValue())
                    .minValue(config.getMinValue())
                    .maxValue(config.getMaxValue())
                    .fileSuffix(config.getFileSuffix())
                    .build())
        .collect(Collectors.toList());
  }

  private List<PipelineOutputDefinition> outputDefinitionsFromConfig(
      List<PipelineOutputDefinitionConfiguration> outputConfigs) {
    if (outputConfigs == null) {
      return java.util.Collections.emptyList();
    }
    return outputConfigs.stream()
        .map(
            config ->
                PipelineOutputDefinition.builder()
                    .name(config.getName())
                    .wdlVariableName(config.getWdlVariableName())
                    .displayName(config.getDisplayName())
                    .description(config.getDescription())
                    .type(config.getType())
                    .isRequired(config.getIsRequired() != null && config.getIsRequired())
                    .build())
        .collect(Collectors.toList());
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

    // TODO ensure that we add a row if one doesn't exist - don't need to getPipeline, can just get
    // pipeline runtime metadata
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

  /**
   * Resolve a pipeline by canonical pipelineKey ({pipeline_name}_v{version}).
   *
   * @param pipelineKey canonical pipeline key
   * @param showHidden whether to include hidden pipelines
   * @return merged pipeline definition/runtime metadata
   */
  public Pipeline getPipelineByPipelineKey(String pipelineKey, boolean showHidden) {
    PipelinesEnum pipelineName = enumFromPipelineKey(pipelineKey);
    Integer pipelineVersion = versionFromPipelineKey(pipelineKey);
    return getPipeline(pipelineName, pipelineVersion, showHidden);
  }
}
