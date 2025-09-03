package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import jakarta.persistence.Column;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.MappedSuperclass;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Getter
@Setter
@NoArgsConstructor
@MappedSuperclass
public abstract class BasePipelineVariableDefinition {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_id")
  private Long pipelineId;

  @Column(name = "name", nullable = false)
  private String name;

  @Column(name = "wdl_variable_name", nullable = false)
  private String wdlVariableName;

  @Column(name = "type", nullable = false)
  private PipelineVariableTypesEnum type;

  @Column(name = "is_required", nullable = false)
  private Boolean isRequired;

  protected BasePipelineVariableDefinition(
      Long pipelineId,
      String name,
      String wdlVariableName,
      PipelineVariableTypesEnum type,
      boolean isRequired) {
    this.pipelineId = pipelineId;
    this.name = name;
    this.wdlVariableName = wdlVariableName;
    this.type = type;
    this.isRequired = isRequired;
  }
}
