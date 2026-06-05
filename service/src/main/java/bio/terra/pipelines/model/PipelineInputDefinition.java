package bio.terra.pipelines.model;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.ToString;
import lombok.experimental.SuperBuilder;

/**
 * Immutable domain model representing a pipeline input definition. Shared variable fields are
 * inherited from BasePipelineVariableDefinition.
 */
@Getter
@ToString(callSuper = true)
@EqualsAndHashCode(callSuper = true)
@NoArgsConstructor(force = true)
@SuperBuilder(toBuilder = true)
public class PipelineInputDefinition extends BasePipelineVariableDefinition {

  private final boolean userProvided;
  private final String defaultValue;
  private final Double minValue;
  private final Double maxValue;
  private final String fileSuffix;

  public static PipelineInputDefinition inputDefinitionFromConfiguration(
      PipelineConfigurations.PipelineInputDefinitionConfiguration config) {
    return PipelineInputDefinition.builder()
        .name(config.getName())
        .wdlVariableName(config.getWdlVariableName())
        .displayName(config.getDisplayName())
        .description(config.getDescription())
        .type(config.getType())
        .isRequired(config.getIsRequired())
        .userProvided(config.getUserProvided())
        .defaultValue(config.getDefaultValue())
        .minValue(config.getMinValue())
        .maxValue(config.getMaxValue())
        .fileSuffix(config.getFileSuffix())
        .build();
  }
}
