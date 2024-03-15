package bio.terra.pipelines.db.entities;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_inputs")
public class PipelineInput {
  @Id
  @Column(name = "job_id", nullable = false, unique = true)
  private Long jobId;

  @Column(name = "inputs", nullable = false)
  private String inputs;
}
