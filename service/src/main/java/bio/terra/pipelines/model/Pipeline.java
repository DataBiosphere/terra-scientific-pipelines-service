package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import java.math.BigDecimal;
import java.time.Instant;
import java.util.Collections;
import java.util.List;
import lombok.Builder;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;

/**
 * Immutable domain model representing a pipeline as seen by the service and controller layers. This
 * combines pipeline definition content (from YAML via PipelineConfigurations) with runtime metadata
 * (from the pipeline_runtime_metadata table).
 */
@Getter
@ToString
@EqualsAndHashCode
@Builder(toBuilder = true)
public class Pipeline {

  // --- Identity (config) ---
  private final PipelinesEnum name;
  private final Integer version;
  private final String pipelineKey;

  // --- Display metadata (from config) ---
  private final String displayName;
  private final String description;
  private final String pipelineType;

  // --- Definition structure (from config) ---
  private final String toolName;
  private final List<PipelineInputDefinition> inputDefinitions;
  private final List<PipelineOutputDefinition> outputDefinitions;

  // --- Execution metadata (from config) ---
  private final BigDecimal memoryRetryMultiplier;

  // --- Runtime metadata (from pipeline_runtime_metadata table) ---
  private final String toolVersion;
  private final String workspaceBillingProject;
  private final String workspaceName;
  private final String workspaceStorageContainerName;
  private final String workspaceGoogleProject;
  private final boolean hidden;
  private final Instant updated;

  /**
   * Returns an immutable copy of the pipeline input definitions.
   *
   * @return immutable list of pipeline input definitions
   */
  public List<PipelineInputDefinition> getInputDefinitions() {
    return Collections.unmodifiableList(inputDefinitions);
  }

  /**
   * Returns an immutable copy of the pipeline output definitions.
   *
   * @return immutable list of pipeline output definitions
   */
  public List<PipelineOutputDefinition> getOutputDefinitions() {
    return Collections.unmodifiableList(outputDefinitions);
  }
}
