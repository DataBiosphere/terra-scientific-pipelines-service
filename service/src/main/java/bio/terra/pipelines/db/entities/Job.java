package bio.terra.pipelines.db.entities;

import java.time.Instant;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Table;

@Entity
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

  public Job() {}

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

  public String getJobId() {
    return jobId;
  }

  public void setJobId(String jobId) {
    this.jobId = jobId;
  }

  public String getUserId() {
    return userId;
  }

  public void setUserId(String userId) {
    this.userId = userId;
  }

  public String getPipelineId() {
    return pipelineId;
  }

  public void setPipelineId(String pipelineId) {
    this.pipelineId = pipelineId;
  }

  public String getPipelineVersion() {
    return pipelineVersion;
  }

  public void setPipelineVersion(String pipelineVersion) {
    this.pipelineVersion = pipelineVersion;
  }

  public Instant getTimeSubmitted() {
    return timeSubmitted;
  }

  public void setTimeSubmitted(Instant timeSubmitted) {
    this.timeSubmitted = timeSubmitted;
  }

  public Instant getTimeCompleted() {
    return timeCompleted;
  }

  public void setTimeCompleted(Instant timeCompleted) {
    this.timeCompleted = timeCompleted;
  }

  public String getStatus() {
    return status;
  }

  public void setStatus(String status) {
    this.status = status;
  }
}
