package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import lombok.Builder;
import lombok.Getter;
import lombok.ToString;

/**
 * Immutable domain model representing a pipeline as seen by the service and controller layers. This
 * combines pipeline definition content (from YAML via PipelineDefinition) with runtime metadata
 * (from the pipeline_runtime_metadata table).
 *
 * <p>This is distinct from the JPA entity {@code bio.terra.pipelines.db.entities.Pipeline}, which
 * is used for DB persistence only. Use this class in service return types, controller logic, and
 * test fixtures that care about pipeline semantics.
 */
@Getter
@ToString
@Builder(toBuilder = true)
public class Pipeline {

  // --- Identity (from PipelineDefinition) ---
  private final PipelinesEnum name;
  private final Integer version;
  private final String pipelineKey;

  // --- Display metadata (from PipelineDefinition) ---
  private final String displayName;
  private final String description;
  private final String pipelineType;

  // --- Definition structure (from PipelineDefinition) ---
  private final String toolName;
  private final List<PipelineInputDefinition> pipelineInputDefinitions;
  private final List<PipelineOutputDefinition> pipelineOutputDefinitions;
  private final PipelineQuota quota;

  // --- Execution metadata (from PipelineDefinition / YAML) ---
  private final List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl;
  private final String storageWorkspaceContainerUrl;
  private final Map<String, String> inputsWithCustomValues;
  private final BigDecimal memoryRetryMultiplier;

  // --- Runtime metadata (from pipeline_runtime_metadata table) ---
  /** Legacy DB id from the pipelines table. May be null if not yet persisted. */
  private final Long id;

  private final String toolVersion;
  private final String workspaceBillingProject;
  private final String workspaceName;
  private final String workspaceStorageContainerName;
  private final String workspaceGoogleProject;

  /**
   * -- GETTER -- Convenience method mirroring the entity Pipeline.isHidden() naming used by many
   * callers.
   *
   * @return whether this pipeline is hidden
   */
  private final boolean hidden;

  /**
   * Returns an immutable copy of the pipeline input definitions.
   *
   * @return immutable list of pipeline input definitions
   */
  public List<PipelineInputDefinition> getPipelineInputDefinitions() {
    return pipelineInputDefinitions != null
        ? Collections.unmodifiableList(pipelineInputDefinitions)
        : Collections.emptyList();
  }

  /**
   * Returns an immutable copy of the pipeline output definitions.
   *
   * @return immutable list of pipeline output definitions
   */
  public List<PipelineOutputDefinition> getPipelineOutputDefinitions() {
    return pipelineOutputDefinitions != null
        ? Collections.unmodifiableList(pipelineOutputDefinitions)
        : Collections.emptyList();
  }

  /**
   * Returns an immutable copy of the input keys to prepend with the storage workspace container
   * URL.
   *
   * @return immutable list of input keys
   */
  public List<String> getInputKeysToPrependWithStorageWorkspaceContainerUrl() {
    return inputKeysToPrependWithStorageWorkspaceContainerUrl != null
        ? Collections.unmodifiableList(inputKeysToPrependWithStorageWorkspaceContainerUrl)
        : Collections.emptyList();
  }

  /**
   * Returns an immutable copy of the inputs with custom values map.
   *
   * @return immutable map of custom value inputs
   */
  public Map<String, String> getInputsWithCustomValues() {
    return inputsWithCustomValues != null
        ? Collections.unmodifiableMap(inputsWithCustomValues)
        : Collections.emptyMap();
  }

  /**
   * Build a Pipeline model from a PipelineDefinition and optional PipelineRuntimeMetadata entity.
   * If runtimeMetadata is null, runtime fields (toolVersion, workspace details) will be null and
   * hidden will default to false.
   *
   * @param definition the pipeline definition from YAML
   * @param runtimeMetadata optional runtime metadata from the database
   * @return the combined Pipeline model
   */
  public static Pipeline fromDefinitionAndRuntime(
      PipelineDefinition definition,
      bio.terra.pipelines.db.entities.PipelineRuntimeMetadata runtimeMetadata) {
    PipelineBuilder builder =
        Pipeline.builder()
            .name(definition.getName())
            .version(definition.getVersion())
            .pipelineKey(definition.getPipelineKey())
            .displayName(definition.getDisplayName())
            .description(definition.getDescription())
            .pipelineType(definition.getPipelineType())
            .toolName(definition.getToolName())
            .pipelineInputDefinitions(definition.getInputs())
            .pipelineOutputDefinitions(definition.getOutputs())
            .quota(definition.getQuota())
            .inputKeysToPrependWithStorageWorkspaceContainerUrl(
                definition.getInputKeysToPrependWithStorageWorkspaceContainerUrl())
            .storageWorkspaceContainerUrl(definition.getStorageWorkspaceContainerUrl())
            .inputsWithCustomValues(definition.getInputsWithCustomValues())
            .memoryRetryMultiplier(definition.getMemoryRetryMultiplier())
            .hidden(false);

    if (runtimeMetadata != null) {
      builder
          .hidden(runtimeMetadata.isHidden())
          .toolVersion(runtimeMetadata.getToolVersion())
          .workspaceBillingProject(runtimeMetadata.getWorkspaceBillingProject())
          .workspaceName(runtimeMetadata.getWorkspaceName())
          .workspaceStorageContainerName(runtimeMetadata.getWorkspaceStorageContainerName())
          .workspaceGoogleProject(runtimeMetadata.getWorkspaceGoogleProject());
    }

    return builder.build();
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    Pipeline that = (Pipeline) o;
    return hidden == that.hidden
        && name == that.name
        && Objects.equals(version, that.version)
        && Objects.equals(pipelineKey, that.pipelineKey)
        && Objects.equals(displayName, that.displayName)
        && Objects.equals(description, that.description)
        && Objects.equals(pipelineType, that.pipelineType)
        && Objects.equals(toolName, that.toolName)
        && Objects.equals(pipelineInputDefinitions, that.pipelineInputDefinitions)
        && Objects.equals(pipelineOutputDefinitions, that.pipelineOutputDefinitions)
        && Objects.equals(quota, that.quota)
        && Objects.equals(
            inputKeysToPrependWithStorageWorkspaceContainerUrl,
            that.inputKeysToPrependWithStorageWorkspaceContainerUrl)
        && Objects.equals(storageWorkspaceContainerUrl, that.storageWorkspaceContainerUrl)
        && Objects.equals(inputsWithCustomValues, that.inputsWithCustomValues)
        && Objects.equals(memoryRetryMultiplier, that.memoryRetryMultiplier)
        && Objects.equals(id, that.id)
        && Objects.equals(toolVersion, that.toolVersion)
        && Objects.equals(workspaceBillingProject, that.workspaceBillingProject)
        && Objects.equals(workspaceName, that.workspaceName)
        && Objects.equals(workspaceStorageContainerName, that.workspaceStorageContainerName)
        && Objects.equals(workspaceGoogleProject, that.workspaceGoogleProject);
  }

  @Override
  public int hashCode() {
    return Objects.hash(
        name,
        version,
        pipelineKey,
        displayName,
        description,
        pipelineType,
        toolName,
        pipelineInputDefinitions,
        pipelineOutputDefinitions,
        quota,
        inputKeysToPrependWithStorageWorkspaceContainerUrl,
        storageWorkspaceContainerUrl,
        inputsWithCustomValues,
        memoryRetryMultiplier,
        id,
        toolVersion,
        workspaceBillingProject,
        workspaceName,
        workspaceStorageContainerName,
        workspaceGoogleProject,
        hidden);
  }
}
