package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import lombok.Builder;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;

/**
 * Immutable domain model representing a pipeline as seen by the service and controller layers. This
 * combines pipeline definition content (from YAML via PipelineDefinition) with runtime metadata
 * (from the pipeline_runtime_metadata table).
 */
@Getter
@ToString
@EqualsAndHashCode
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
  private final List<PipelineInputDefinition> inputDefinitions;
  private final List<PipelineOutputDefinition> outputDefinitions;
  private final PipelineQuota quota;

  // --- Execution metadata (from PipelineDefinition / YAML) ---
  private final List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl;
  private final String storageWorkspaceContainerUrl;
  private final Map<String, String> inputsWithCustomValues;
  private final BigDecimal memoryRetryMultiplier;

  // --- Runtime metadata (from pipeline_runtime_metadata table) ---
  private final String toolVersion;
  private final String workspaceBillingProject;
  private final String workspaceName;
  private final String workspaceStorageContainerName;
  private final String workspaceGoogleProject;
  private final boolean hidden;

  /**
   * Returns an immutable copy of the pipeline input definitions.
   *
   * @return immutable list of pipeline input definitions
   */
  public List<PipelineInputDefinition> getInputDefinitions() {
    return unmodifiableOrEmpty(inputDefinitions);
  }

  /**
   * Returns an immutable copy of the pipeline output definitions.
   *
   * @return immutable list of pipeline output definitions
   */
  public List<PipelineOutputDefinition> getOutputDefinitions() {
    return unmodifiableOrEmpty(outputDefinitions);
  }

  /**
   * Returns an immutable copy of the input keys to prepend with the storage workspace container
   * URL.
   *
   * @return immutable list of input keys
   */
  public List<String> getInputKeysToPrependWithStorageWorkspaceContainerUrl() {
    return unmodifiableOrEmpty(inputKeysToPrependWithStorageWorkspaceContainerUrl);
  }

  /**
   * Returns an immutable copy of the inputs with custom values map.
   *
   * @return immutable map of custom value inputs
   */
  public Map<String, String> getInputsWithCustomValues() {
    return unmodifiableOrEmpty(inputsWithCustomValues);
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
      PipelineDefinition definition, PipelineRuntimeMetadata runtimeMetadata) {
    PipelineBuilder builder =
        Pipeline.builder()
            .name(definition.getName())
            .version(definition.getVersion())
            .pipelineKey(definition.getPipelineKey())
            .displayName(definition.getDisplayName())
            .description(definition.getDescription())
            .pipelineType(definition.getPipelineType())
            .toolName(definition.getToolName())
            .inputDefinitions(definition.getInputDefinitions())
            .outputDefinitions(definition.getOutputDefinitions())
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

  private static <T> List<T> unmodifiableOrEmpty(List<T> values) {
    return values != null ? Collections.unmodifiableList(values) : Collections.emptyList();
  }

  private static <K, V> Map<K, V> unmodifiableOrEmpty(Map<K, V> values) {
    return values != null ? Collections.unmodifiableMap(values) : Collections.emptyMap();
  }
}
