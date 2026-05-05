package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Getter
@Setter
@NoArgsConstructor
public class PipelineOutputDefinition extends PipelineVariableDefinition {
  public PipelineOutputDefinition(
      String pipelineKey,
      String name,
      String wdlVariableName,
      String displayName,
      String description,
      PipelineVariableTypesEnum type,
      boolean isRequired) {
    super(pipelineKey, name, wdlVariableName, displayName, description, type, isRequired);
  }
}
