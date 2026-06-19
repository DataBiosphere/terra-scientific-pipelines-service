package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import lombok.Builder;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;

/** Shared config and domain model representing a pipeline input definition. */
@Getter
@ToString
@EqualsAndHashCode
@Builder
public class PipelineInputDefinition {
  private final String name;
  private final String wdlVariableName;
  private final String displayName;
  private final String description;
  private final PipelineVariableTypesEnum type;
  private final boolean isRequired;
  private final boolean userProvided;
  private final String defaultValue;
  private final Double minValue;
  private final Double maxValue;
  private final String fileSuffix;
  private final String validationRegex;
}
