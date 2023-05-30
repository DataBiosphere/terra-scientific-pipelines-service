package bio.terra.pipelines.service.model;

import bio.terra.pipelines.db.DbJob;
import java.time.Instant;
import java.util.Optional;
import java.util.UUID;

public record Job(
    UUID jobId,
    String userId,
    String pipelineId,
    String pipelineVersion,
    Instant timeSubmitted,
    Optional<Instant> timeCompleted,
    String status) {

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
