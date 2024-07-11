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
@Table(name = "pipeline_outputs")
public class PipelineOutputs {
  @Id
  @Column(name = "job_id", nullable = false, unique = true)
  private Long jobId;

  @Column(name = "outputs", nullable = false)
  private String outputs;
}
