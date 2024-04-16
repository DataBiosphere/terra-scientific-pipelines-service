package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_inputs_definitions")
public class PipelineInputsDefinition {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_id")
  private Long pipelineId;

  @Column(name = "input_name", nullable = false)
  private String inputName;

  @Column(name = "input_type", nullable = false)
  private String inputType;

  @Column(name = "is_required", nullable = false)
  private Boolean isRequired;

  public PipelineInputsDefinition(
      Long pipelineId, String inputName, String inputType, Boolean isRequired) {
    this.pipelineId = pipelineId;
    this.inputName = inputName;
    this.inputType = inputType;
    this.isRequired = isRequired;
  }
}
