package bio.terra.pipelines.db.entities;

import java.util.UUID;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Table;
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
  private UUID jobId;

  @Column(name = "inputs", nullable = false)
  private String inputs;
}
