package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import java.time.LocalDateTime;
import java.util.UUID;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.hibernate.annotations.CreationTimestamp;
import org.hibernate.annotations.SourceType;
import org.hibernate.annotations.UpdateTimestamp;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_runs")
@SuppressWarnings("java:S107") // Disable "Methods should not have too many parameters"
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
  @CreationTimestamp(source = SourceType.DB)
  private LocalDateTime created;

  @Column(name = "updated", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private LocalDateTime updated;

  @Column(name = "status", nullable = false)
  private String status;

  @Column(name = "description")
  private String description;

  @Column(name = "result_url")
  private String resultUrl;

  @Column(name = "is_success")
  private Boolean isSuccess;

  @Column(name = "output")
  private String output;

  /** Constructor for in progress or complete PipelineRun. */
  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      LocalDateTime created,
      LocalDateTime updated,
      String status,
      String description,
      String resultUrl,
      Boolean isSuccess,
      String output) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.created = created;
    this.updated = updated;
    this.status = status;
    this.description = description;
    this.resultUrl = resultUrl;
    this.isSuccess = isSuccess;
    this.output = output;
  }

  /** Constructor for creating a new pipeline run. Timestamps are auto-generated. */
  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      String status,
      String description,
      String resultUrl) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.status = status;
    this.description = description;
    this.resultUrl = resultUrl;
  }
}