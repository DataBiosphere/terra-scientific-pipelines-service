package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
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
@Table(name = "pipeline_input_definitions")
public class PipelineInputDefinition extends BasePipelineVariableDefinition {
  @Column(name = "is_required", nullable = false)
  private Boolean isRequired;

  @Column(name = "user_provided", nullable = false)
  private Boolean userProvided;

  @Column(name = "default_value")
  private String defaultValue; // must be a String representation of the value

  @Column(name = "type", nullable = false)
  private PipelineInputTypesEnum type;

  public PipelineInputDefinition(
      Long pipelineId,
      String name,
      String wdlVariableName,
      PipelineInputTypesEnum type,
      Boolean isRequired,
      Boolean userProvided,
      String defaultValue) {
    super(pipelineId, name, wdlVariableName);
    this.type = type;
    this.isRequired = isRequired;
    this.userProvided = userProvided;
    this.defaultValue = defaultValue;
  }
}
