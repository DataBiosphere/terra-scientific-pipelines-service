package bio.terra.pipelines.db.entities;

import java.time.Instant;
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
@Table(name = "jobs")
public class Job {
  @Id
  @Column(name = "job_id", nullable = false)
  private String jobId;

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
      String jobId,
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
