package bio.terra.pipelines.db.entities;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_output_definitions")
public class PipelineOutputDefinition {
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
  private String type;

  public PipelineOutputDefinition(
      Long pipelineId, String name, String wdlVariableName, String type) {
    this.pipelineId = pipelineId;
    this.name = name;
    this.wdlVariableName = wdlVariableName;
    this.type = type;
  }
}
