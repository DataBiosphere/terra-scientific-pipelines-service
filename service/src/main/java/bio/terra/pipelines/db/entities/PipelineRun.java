package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import java.time.LocalDateTime;
import java.util.UUID;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_runs")
public class PipelineRun {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "job_id", nullable = false, unique = true)
  private UUID jobId;

  @Column(name = "user_id", nullable = false)
  private String userId;

  @Column(name = "pipeline_id", nullable = false)
  private Long pipelineId;

  @Column(name = "created", insertable = false)
  private LocalDateTime created;

  @Column(name = "updated", insertable = false)
  private LocalDateTime updated;

  @Column(name = "status", nullable = false)
  private String status;

  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      LocalDateTime created,
      LocalDateTime updated,
      String status) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.created = created;
    this.updated = updated;
    this.status = status;
  }

  public PipelineRun(UUID jobId, String userId, Long pipelineId, String status) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.status = status;
  }
}
