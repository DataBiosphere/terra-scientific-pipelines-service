package bio.terra.pipelines.model;

import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;
import lombok.experimental.SuperBuilder;

/**
 * Immutable domain model representing a pipeline output definition. Shared variable fields are
 * inherited from BasePipelineVariableDefinition.
 */
@Getter
@ToString(callSuper = true)
@EqualsAndHashCode(callSuper = true)
@SuperBuilder(toBuilder = true)
public class PipelineOutputDefinition extends BasePipelineVariableDefinition {}
