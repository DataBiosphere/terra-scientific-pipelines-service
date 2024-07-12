package bio.terra.pipelines.db.entities;

import jakarta.persistence.Column;
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

  @Column(name = "type", nullable = false)
  private String type;

  public PipelineOutputDefinition(
      Long pipelineId, String name, String wdlVariableName, String type) {
    super(pipelineId, name, wdlVariableName);
    this.type = type;
  }
}
