package bio.terra.pipelines.service.model;

import bio.terra.pipelines.db.DbJob;
import java.time.Instant;
import java.util.Optional;
import java.util.UUID;

public class Job {
  private final UUID jobId;
  private final String userId;
  private final String pipelineId;
  private final String pipelineVersion;
  private final Instant timeSubmitted;
  private final Optional<Instant> timeCompleted;
  private final String status;
  private final Object pipelineInputs;

  public Job(
      UUID jobId,
      String userId,
      String pipelineId,
      String pipelineVersion,
      Instant timeSubmitted,
      Optional<Instant> timeCompleted,
      String status,
      Object pipelineInputs) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.pipelineVersion = pipelineVersion;
    this.timeSubmitted = timeSubmitted;
    this.timeCompleted = timeCompleted;
    this.status = status;
    this.pipelineInputs = pipelineInputs;
  }

  public UUID getJobId() {
    return jobId;
  }

  public String getUserId() {
    return userId;
  }

  public String getPipelineId() {
    return pipelineId;
  }

  public String getPipelineVersion() {
    return pipelineVersion;
  }

  public Instant getTimeSubmitted() {
    return timeSubmitted;
  }

  public Optional<Instant> getTimeCompleted() {
    return timeCompleted;
  }

  public String getStatus() {
    return status;
  }

  public Object getPipelineInputs() {
    return pipelineInputs;
  }

  public static Job fromDb(DbJob dbJob) {
    return new Job(
        dbJob.jobId(),
        dbJob.userId(),
        dbJob.pipelineId(),
        dbJob.pipelineVersion(),
        dbJob.timeSubmitted(),
        dbJob.timeCompleted(),
        dbJob.status(),
        dbJob.pipelineInputs());
  }
}
