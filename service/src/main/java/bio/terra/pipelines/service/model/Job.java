package bio.terra.pipelines.service.model;

import bio.terra.pipelines.db.entities.DbJob;
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

  public Job(
      UUID jobId,
      String userId,
      String pipelineId,
      String pipelineVersion,
      Instant timeSubmitted,
      Optional<Instant> timeCompleted,
      String status) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.pipelineVersion = pipelineVersion;
    this.timeSubmitted = timeSubmitted;
    this.timeCompleted = timeCompleted;
    this.status = status;
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

  public static Job fromDb(DbJob dbJob) {
    return new Job(
        dbJob.getJobId(),
        dbJob.getUserId(),
        dbJob.getPipelineId(),
        dbJob.getPipelineVersion(),
        dbJob.getTimeSubmitted(),
        Optional.of(dbJob.getTimeCompleted()),
        dbJob.getStatus());
  }
}
