package bio.terra.pipelines.model;

import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;
import lombok.experimental.SuperBuilder;

/**
 * Immutable domain model representing a pipeline input definition. Shared variable fields are
 * inherited from BasePipelineVariableDefinition.
 */
@Getter
@ToString(callSuper = true)
@EqualsAndHashCode(callSuper = true)
@SuperBuilder(toBuilder = true)
public class PipelineInputDefinition extends BasePipelineVariableDefinition {

  private final boolean userProvided;
  private final boolean expectsCustomValue;
  private final String defaultValue;
  private final Double minValue;
  private final Double maxValue;
  private final String fileSuffix;
}
