package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import lombok.AllArgsConstructor;
import lombok.Builder;
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
  private final List<PipelineInputDefinition> inputs;
  private final List<PipelineOutputDefinition> outputs;
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
  public List<PipelineInputDefinition> getInputs() {
    return inputs != null ? Collections.unmodifiableList(inputs) : Collections.emptyList();
  }

  /**
   * Returns an immutable copy of the outputs list.
   *
   * @return immutable list of pipeline outputs
   */
  public List<PipelineOutputDefinition> getOutputs() {
    return outputs != null ? Collections.unmodifiableList(outputs) : Collections.emptyList();
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

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    PipelineDefinition that = (PipelineDefinition) o;
    return name == that.name
        && Objects.equals(version, that.version)
        && Objects.equals(pipelineKey, that.pipelineKey)
        && Objects.equals(displayName, that.displayName)
        && Objects.equals(description, that.description)
        && Objects.equals(pipelineType, that.pipelineType)
        && Objects.equals(toolName, that.toolName)
        && Objects.equals(inputs, that.inputs)
        && Objects.equals(outputs, that.outputs)
        && Objects.equals(quota, that.quota)
        && Objects.equals(
            inputKeysToPrependWithStorageWorkspaceContainerUrl,
            that.inputKeysToPrependWithStorageWorkspaceContainerUrl)
        && Objects.equals(storageWorkspaceContainerUrl, that.storageWorkspaceContainerUrl)
        && Objects.equals(inputsWithCustomValues, that.inputsWithCustomValues)
        && Objects.equals(memoryRetryMultiplier, that.memoryRetryMultiplier);
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
        inputs,
        outputs,
        quota,
        inputKeysToPrependWithStorageWorkspaceContainerUrl,
        storageWorkspaceContainerUrl,
        inputsWithCustomValues,
        memoryRetryMultiplier);
  }
}
