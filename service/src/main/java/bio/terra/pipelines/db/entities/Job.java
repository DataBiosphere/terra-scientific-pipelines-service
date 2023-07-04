package bio.terra.pipelines.db.entities;

import java.time.Instant;
import java.util.UUID;
import javax.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "jobs")
public class Job {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "job_id", nullable = false, unique = true)
  private UUID jobId;

  @Column(name = "user_id", nullable = false)
  private String userId;

  @Column(name = "pipeline_id", nullable = false)
  private String pipelineId;

  @Column(name = "pipeline_version", nullable = false)
  private String pipelineVersion;

  @Column(name = "time_submitted", nullable = false)
  private Instant timeSubmitted;

  @Column(name = "time_completed")
  private Instant timeCompleted;

  @Column(name = "status", nullable = false)
  private String status;

  public Job(
      UUID jobId,
      String userId,
      String pipelineId,
      String pipelineVersion,
      Instant timeSubmitted,
      Instant timeCompleted,
      String status) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.pipelineVersion = pipelineVersion;
    this.timeSubmitted = timeSubmitted;
    this.timeCompleted = timeCompleted;
    this.status = status;
  }
}
