package bio.terra.pipelines.db.entities;

import java.util.UUID;
import javax.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "imputation_jobs")
public class ImputationJob {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "job_id", nullable = false, unique = true)
  private UUID jobId;

  @Column(name = "user_id", nullable = false)
  private String userId;

  @Column(name = "pipeline_version", nullable = false)
  private String pipelineVersion;

  public ImputationJob(UUID jobId, String userId, String pipelineVersion) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineVersion = pipelineVersion;
  }
}