package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.ToString;
import lombok.experimental.SuperBuilder;

/** Shared immutable fields for pipeline input/output definitions in the domain model. */
@Getter
@ToString
@EqualsAndHashCode
@NoArgsConstructor(force = true)
@SuperBuilder(toBuilder = true)
public abstract class BasePipelineVariableDefinition {
  private final String name;
  private final String wdlVariableName;
  private final String displayName;
  private final String description;
  private final PipelineVariableTypesEnum type;
  private final boolean isRequired;
}
