package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
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
  private final PipelinesRepository pipelinesRepository;
  private final RawlsService rawlsService;
  private final SamService samService;

  // \\d means digit [0-9]; \\w means word-like [a-zA-Z0-9_]
  private static final String SEM_VER_REGEX_STRING = "^(\\d+|\\w*v\\d+)\\.(\\d+)\\.(\\d+)$";

  @Autowired
  public PipelinesService(
      PipelineDefinitionProvider pipelineDefinitionProvider,
      PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository,
      PipelinesRepository pipelinesRepository,
      RawlsService rawlsService,
      SamService samService) {
    this.pipelineDefinitionProvider = pipelineDefinitionProvider;
    this.pipelineRuntimeMetadataRepository = pipelineRuntimeMetadataRepository;
    this.pipelinesRepository = pipelinesRepository;
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
        // Try to load runtime metadata from database
        PipelineRuntimeMetadata runtimeMetadata =
            pipelineRuntimeMetadataRepository.findById(definition.getPipelineKey()).orElse(null);

        // Apply hidden filter
        if (runtimeMetadata != null && runtimeMetadata.isHidden() && !showHidden) {
          continue;
        }

        // Construct Pipeline entity from definition and runtime metadata
        Pipeline pipeline = constructPipelineFromDefinition(definition, runtimeMetadata);
        pipelines.add(pipeline);
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

    // Load definition from YAML
    PipelineDefinition definition =
        pipelineDefinitionProvider.getPipelineDefinition(pipelineName, pipelineVersion);

    // Try to load runtime metadata from database
    PipelineRuntimeMetadata runtimeMetadata =
        pipelineRuntimeMetadataRepository.findById(definition.getPipelineKey()).orElse(null);

    // Check hidden status
    if (runtimeMetadata != null && runtimeMetadata.isHidden() && !showHidden) {
      throw new NotFoundException(
          "Pipeline not found for pipelineName %s and version %s"
              .formatted(pipelineName, pipelineVersion));
    }

    return constructPipelineFromDefinition(definition, runtimeMetadata);
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

    // Find first non-hidden version if showHidden is false
    for (PipelineDefinition definition : definitions) {
      PipelineRuntimeMetadata runtimeMetadata =
          pipelineRuntimeMetadataRepository.findById(definition.getPipelineKey()).orElse(null);

      if (runtimeMetadata != null && runtimeMetadata.isHidden() && !showHidden) {
        continue;
      }

      return constructPipelineFromDefinition(definition, runtimeMetadata);
    }

    throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
  }

  /**
   * Get a pipeline by its database ID. This method is kept for backward compatibility with existing
   * code paths that may reference pipelines by ID.
   *
   * @param pipelineId the pipeline ID
   * @return the pipeline
   * @throws NotFoundException if no pipeline exists with that ID
   * @deprecated only used for tests, should be removed once tests are refactored
   */
  @Deprecated(since = "2.8.0")
  public Pipeline getPipelineById(Long pipelineId) {
    logger.info("Get a specific pipeline for pipelineId {}", pipelineId);
    Pipeline dbResult = pipelinesRepository.findById(pipelineId).orElse(null);
    if (dbResult == null) {
      throw new NotFoundException("Pipeline not found for pipelineId %s".formatted(pipelineId));
    }
    if (dbResult.getPipelineKey() == null || dbResult.getPipelineKey().isBlank()) {
      return dbResult;
    }
    return getPipelineByKey(dbResult.getPipelineKey());
  }

  /**
   * Resolve a pipeline by canonical pipelineKey ({pipeline_name}_v{version}).
   *
   * @param pipelineKey canonical pipeline key
   * @return merged pipeline definition/runtime metadata
   */
  public Pipeline getPipelineByKey(String pipelineKey) {
    int separatorIndex = pipelineKey.lastIndexOf("_v");
    if (separatorIndex < 1 || separatorIndex >= pipelineKey.length() - 2) {
      throw new NotFoundException("Pipeline not found for pipelineKey %s".formatted(pipelineKey));
    }

    String pipelineNameValue = pipelineKey.substring(0, separatorIndex);
    Integer pipelineVersion;
    try {
      pipelineVersion = Integer.valueOf(pipelineKey.substring(separatorIndex + 2));
    } catch (NumberFormatException e) {
      throw new NotFoundException("Pipeline not found for pipelineKey %s".formatted(pipelineKey));
    }

    PipelinesEnum pipelineName =
        java.util.Arrays.stream(PipelinesEnum.values())
            .filter(p -> p.getValue().equals(pipelineNameValue))
            .findFirst()
            .orElseThrow(
                () ->
                    new NotFoundException(
                        "Pipeline not found for pipelineKey %s".formatted(pipelineKey)));

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

  /**
   * Construct a Pipeline entity from a PipelineDefinition and optional runtime metadata.
   *
   * <p>This bridges between the YAML-backed domain model (PipelineDefinition) and the database
   * entity (Pipeline). Definition data comes from YAML; runtime metadata comes from the database
   * (or defaults if not yet persisted).
   *
   * @param definition the pipeline definition from YAML
   * @param runtimeMetadata optional runtime metadata from database (may be null)
   * @return a Pipeline entity with merged definition and runtime data
   */
  private Pipeline constructPipelineFromDefinition(
      PipelineDefinition definition, PipelineRuntimeMetadata runtimeMetadata) {
    Pipeline pipeline = new Pipeline();

    // Identity
    pipeline.setName(definition.getName());
    pipeline.setVersion(definition.getVersion());
    pipeline.setPipelineKey(definition.getPipelineKey());
    // Compatibility-only: pipeline_id still exists in pipeline_runs and Stairway map.
    Pipeline legacyPipelineRow =
        pipelinesRepository.findByNameAndVersion(definition.getName(), definition.getVersion());
    if (legacyPipelineRow != null) {
      pipeline.setId(legacyPipelineRow.getId());
    }

    // Definition metadata
    pipeline.setDisplayName(definition.getDisplayName());
    pipeline.setDescription(definition.getDescription());
    pipeline.setPipelineType(definition.getPipelineType());
    pipeline.setToolName(definition.getToolName());

    // Input/output definitions from YAML
    List<PipelineInputDefinition> inputDefs = new ArrayList<>();
    for (var domainInput : definition.getInputs()) {
      PipelineInputDefinition dbInput = new PipelineInputDefinition();
      dbInput.setName(domainInput.getName());
      dbInput.setWdlVariableName(domainInput.getWdlVariableName());
      dbInput.setDisplayName(domainInput.getDisplayName());
      dbInput.setDescription(domainInput.getDescription());
      dbInput.setType(domainInput.getType());
      dbInput.setRequired(domainInput.isRequired());
      dbInput.setUserProvided(domainInput.isUserProvided());
      dbInput.setExpectsCustomValue(domainInput.isExpectsCustomValue());
      dbInput.setDefaultValue(domainInput.getDefaultValue());
      dbInput.setMinValue(domainInput.getMinValue());
      dbInput.setMaxValue(domainInput.getMaxValue());
      dbInput.setFileSuffix(domainInput.getFileSuffix());
      inputDefs.add(dbInput);
    }
    pipeline.setPipelineInputDefinitions(inputDefs);

    List<PipelineOutputDefinition> outputDefs = new ArrayList<>();
    for (var domainOutput : definition.getOutputs()) {
      PipelineOutputDefinition dbOutput = new PipelineOutputDefinition();
      dbOutput.setName(domainOutput.getName());
      dbOutput.setWdlVariableName(domainOutput.getWdlVariableName());
      dbOutput.setDisplayName(domainOutput.getDisplayName());
      dbOutput.setDescription(domainOutput.getDescription());
      dbOutput.setType(domainOutput.getType());
      dbOutput.setRequired(domainOutput.isRequired());
      outputDefs.add(dbOutput);
    }
    pipeline.setPipelineOutputDefinitions(outputDefs);

    // Runtime metadata (from database, or defaults)
    if (runtimeMetadata != null) {
      pipeline.setHidden(runtimeMetadata.isHidden());
      pipeline.setToolVersion(runtimeMetadata.getToolVersion());
      pipeline.setWorkspaceBillingProject(runtimeMetadata.getWorkspaceBillingProject());
      pipeline.setWorkspaceName(runtimeMetadata.getWorkspaceName());
      pipeline.setWorkspaceStorageContainerName(runtimeMetadata.getWorkspaceStorageContainerName());
      pipeline.setWorkspaceGoogleProject(runtimeMetadata.getWorkspaceGoogleProject());
    } else {
      // Default values for runtime metadata if not persisted yet
      pipeline.setHidden(false);
    }

    return pipeline;
  }

  /**
   * This method is meant to only be called by an admin endpoint to update pipeline parameters such
   * as control workspace information and wdl method version.
   *
   * <p>Calls Rawls to fetch control workspace metadata based on the workspaceBillingProject and
   * workspaceName.
   *
   * <p>Note: definition data (displayName, description, inputs, outputs) is read from YAML and
   * cannot be updated through this method. This method only updates runtime metadata (toolVersion,
   * workspace details, hidden status).
   *
   * @param pipelineName - name of pipeline to update
   * @param pipelineVersion - version of pipeline to update
   * @param isHidden - whether the pipeline should be hidden, if null, will not update the
   *     pipeline's value
   * @param workspaceBillingProject - workspace billing project to update to
   * @param workspaceName - workspace name to update to
   * @param toolVersion - version of the tool expected to run for corresponding pipeline. must align
   *     with pipeline version
   */
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

    // Load the full pipeline (definition + runtime metadata)
    Pipeline pipeline = getPipeline(pipelineName, pipelineVersion, true);

    // Update runtime metadata only
    pipeline.setWorkspaceBillingProject(workspaceBillingProject);
    pipeline.setWorkspaceName(workspaceName);
    pipeline.setWorkspaceStorageContainerName(workspaceStorageContainerUrl);
    pipeline.setWorkspaceGoogleProject(workspaceGoogleProject);
    if (isHidden != null) {
      pipeline.setHidden(isHidden);
    }

    // ensure toolVersion follows semantic versioning regex (can be preceded by a string ending
    // in v)
    final Pattern pattern = Pattern.compile(SEM_VER_REGEX_STRING);
    final Matcher matcher = pattern.matcher(toolVersion);
    if (!matcher.matches()) {
      throw new ValidationException(
          String.format(
              "toolVersion %s does not follow semantic versioning regex %s",
              toolVersion, SEM_VER_REGEX_STRING));
    }
    pipeline.setToolVersion(toolVersion);

    PipelineRuntimeMetadata runtimeMetadata = new PipelineRuntimeMetadata();
    runtimeMetadata.setPipelineKey(pipeline.getPipelineKey());
    runtimeMetadata.setHidden(pipeline.isHidden());
    runtimeMetadata.setToolVersion(pipeline.getToolVersion());
    runtimeMetadata.setWorkspaceBillingProject(pipeline.getWorkspaceBillingProject());
    runtimeMetadata.setWorkspaceName(pipeline.getWorkspaceName());
    runtimeMetadata.setWorkspaceStorageContainerName(pipeline.getWorkspaceStorageContainerName());
    runtimeMetadata.setWorkspaceGoogleProject(pipeline.getWorkspaceGoogleProject());
    pipelineRuntimeMetadataRepository.save(runtimeMetadata);

    // Compatibility mirror while legacy code paths still depend on pipelines table runtime fields.
    Pipeline legacyPipelineRow =
        pipelinesRepository.findByNameAndVersion(pipelineName, pipelineVersion);
    if (legacyPipelineRow != null) {
      legacyPipelineRow.setHidden(pipeline.isHidden());
      legacyPipelineRow.setToolVersion(pipeline.getToolVersion());
      legacyPipelineRow.setWorkspaceBillingProject(pipeline.getWorkspaceBillingProject());
      legacyPipelineRow.setWorkspaceName(pipeline.getWorkspaceName());
      legacyPipelineRow.setWorkspaceStorageContainerName(
          pipeline.getWorkspaceStorageContainerName());
      legacyPipelineRow.setWorkspaceGoogleProject(pipeline.getWorkspaceGoogleProject());
      pipelinesRepository.save(legacyPipelineRow);
    }

    return pipeline;
  }
}
