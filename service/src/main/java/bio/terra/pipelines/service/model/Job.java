package bio.terra.pipelines.service.model;

import bio.terra.pipelines.db.DbJob;
import java.sql.Timestamp;
import java.util.Optional;
import java.util.UUID;

public class Job {
  private final UUID jobId;
  private final String userId;
  private final String pipelineId;
  private final String pipelineVersion;
  private final Timestamp timeSubmitted;
  private final Optional<Timestamp> timeCompleted;
  private final String status;

  public Job(
      UUID jobId,
      String userId,
      String pipelineId,
      String pipelineVersion,
      Timestamp timeSubmitted,
      Optional<Timestamp> timeCompleted,
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

  public Timestamp getTimeSubmitted() {
    return timeSubmitted;
  }

  public Optional<Timestamp> getTimeCompleted() {
    return timeCompleted;
  }

  public String getStatus() {
    return status;
  }

  public static Job fromDb(DbJob dbJob) {
    return new Job(
        dbJob.jobId(),
        dbJob.userId(),
        dbJob.pipelineId(),
        dbJob.pipelineVersion(),
        dbJob.timeSubmitted(),
        dbJob.timeCompleted(),
        dbJob.status());
  }
}
