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
    return new Job.Builder()
        .setJobId(dbJob.jobId())
        .setUserId(dbJob.userId())
        .setPipelineId(dbJob.pipelineId())
        .setPipelineVersion(dbJob.pipelineVersion())
        .setTimeSubmitted(dbJob.timeSubmitted())
        .setTimeCompleted(Optional.of(dbJob.timeCompleted()))
        .setStatus(dbJob.status())
        .build();
  }

  public static class Builder {
    private UUID jobId;
    private String userId;
    private String pipelineId;
    private String pipelineVersion;
    private Timestamp timeSubmitted;
    private Optional<Timestamp> timeCompleted;
    private String status;

    public Builder setJobId(UUID JobId) {
      this.jobId = JobId;
      return this;
    }

    public Builder setUserId(String userId) {
      this.userId = userId;
      return this;
    }

    public Builder setPipelineId(String pipelineId) {
      this.pipelineId = pipelineId;
      return this;
    }

    public Builder setPipelineVersion(String pipelineVersion) {
      this.pipelineVersion = pipelineVersion;
      return this;
    }

    public Builder setTimeSubmitted(Timestamp timeSubmitted) {
      this.timeSubmitted = timeSubmitted;
      return this;
    }

    public Builder setTimeCompleted(Optional<Timestamp> timeCompleted) {
      this.timeCompleted = timeCompleted;
      return this;
    }

    public Builder setStatus(String status) {
      this.status = status;
      return this;
    }

    public Job build() {
      return new Job(
          jobId, userId, pipelineId, pipelineVersion, timeSubmitted, timeCompleted, status);
    }
  }
}
