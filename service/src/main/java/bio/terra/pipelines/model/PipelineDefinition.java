package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;

/**
 * Immutable domain model representing a complete pipeline definition. This is the semantic
 * representation of "what is this pipeline?" and is derived from the YAML pipeline configuration.
 *
 * <p>This is distinct from runtime metadata (like workspace details or toolVersion) which is
 * managed in the database. This class captures the structural and semantic aspects of the pipeline
 * definition.
 */
@Getter
@ToString
@EqualsAndHashCode
@AllArgsConstructor
@Builder(toBuilder = true)
public class PipelineDefinition {

  // Structural identity
  private final PipelinesEnum name;
  private final Integer version;
  private final String pipelineKey; // Format: {pipelineName}_v{version}

  // Display metadata
  private final String displayName;
  private final String description;
  private final String pipelineType;

  // Definition structure
  private final String toolName;
  private final List<PipelineInputDefinition> inputDefinitions;
  private final List<PipelineOutputDefinition> outputDefinitions;
  private final PipelineQuota quota;

  // Execution metadata from YAML
  private final List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl;
  private final String storageWorkspaceContainerUrl;
  private final Map<String, String> inputsWithCustomValues;
  private final BigDecimal memoryRetryMultiplier;

  /**
   * Returns an immutable copy of the inputs list.
   *
   * @return immutable list of pipeline inputs
   */
  public List<PipelineInputDefinition> getInputDefinitions() {
    return inputDefinitions != null
        ? Collections.unmodifiableList(inputDefinitions)
        : Collections.emptyList();
  }

  /**
   * Returns an immutable copy of the outputs list.
   *
   * @return immutable list of pipeline outputs
   */
  public List<PipelineOutputDefinition> getOutputDefinitions() {
    return outputDefinitions != null
        ? Collections.unmodifiableList(outputDefinitions)
        : Collections.emptyList();
  }

  /**
   * Returns an immutable copy of the inputsWithCustomValues map.
   *
   * @return immutable map of custom value inputs
   */
  public Map<String, String> getInputsWithCustomValues() {
    return inputsWithCustomValues != null
        ? Collections.unmodifiableMap(inputsWithCustomValues)
        : Collections.emptyMap();
  }

  /**
   * Returns an immutable copy of the inputKeysToPrependWithStorageWorkspaceContainerUrl list.
   *
   * @return immutable list of input keys to prepend with storage URL
   */
  public List<String> getInputKeysToPrependWithStorageWorkspaceContainerUrl() {
    return inputKeysToPrependWithStorageWorkspaceContainerUrl != null
        ? Collections.unmodifiableList(inputKeysToPrependWithStorageWorkspaceContainerUrl)
        : Collections.emptyList();
  }
}
