package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import jakarta.persistence.Entity;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_output_definitions")
public class PipelineOutputDefinition extends BasePipelineVariableDefinition {

  public PipelineOutputDefinition(
      Long pipelineId, String name, String wdlVariableName, PipelineVariableTypesEnum type) {
    super(pipelineId, name, wdlVariableName, type);
  }
}
