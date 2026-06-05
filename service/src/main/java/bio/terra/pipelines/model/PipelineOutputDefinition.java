package bio.terra.pipelines.model;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.ToString;
import lombok.experimental.SuperBuilder;

/**
 * Immutable domain model representing a pipeline output definition. Shared variable fields are
 * inherited from BasePipelineVariableDefinition.
 */
@Getter
@ToString(callSuper = true)
@EqualsAndHashCode(callSuper = true)
@NoArgsConstructor(force = true)
@SuperBuilder(toBuilder = true)
public class PipelineOutputDefinition extends BasePipelineVariableDefinition {

  public static PipelineOutputDefinition outputDefinitionFromConfiguration(
      PipelineConfigurations.PipelineOutputDefinitionConfiguration config) {
    return PipelineOutputDefinition.builder()
        .name(config.getName())
        .wdlVariableName(config.getWdlVariableName())
        .displayName(config.getDisplayName())
        .description(config.getDescription())
        .type(config.getType())
        .isRequired(config.getIsRequired())
        .build();
  }
}
