package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_input_definitions")
public class PipelineInputDefinition {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_id")
  private Long pipelineId;

  @Column(name = "name", nullable = false)
  private String name;

  @Column(name = "type", nullable = false)
  private String type;

  @Column(name = "is_required", nullable = false)
  private Boolean isRequired;

  public PipelineInputDefinition(Long pipelineId, String name, String type, Boolean isRequired) {
    this.pipelineId = pipelineId;
    this.name = name;
    this.type = type;
    this.isRequired = isRequired;
  }
}
