package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(
    name = "pipeline_outputs",
    uniqueConstraints = {
      @UniqueConstraint(
          name = "pipeline_runs_id_output_name_uk",
          columnNames = {"pipeline_runs_id", "output_name"})
    },
    indexes = {
      @Index(name = "pipeline_outputs_pipeline_runs_id_idx", columnList = "pipeline_runs_id")
    })
public class PipelineOutput {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_runs_id", nullable = false)
  private Long pipelineRunsId;

  @Column(name = "output_name", nullable = false, columnDefinition = "TEXT")
  private String outputName;

  @Column(name = "output_value", nullable = false, columnDefinition = "TEXT")
  private String outputValue;
}
