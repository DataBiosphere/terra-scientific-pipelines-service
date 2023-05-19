package bio.terra.pipelines.db;

import java.time.Instant;
import java.util.Optional;
import java.util.UUID;

/** Record to hold a Job record when processing in the JobDao */
public record DbJob(
    UUID jobId,
    String userId,
    String pipelineId,
    String pipelineVersion,
    Instant timeSubmitted,
    Optional<Instant> timeCompleted,
    String status,
    Object pipelineInputs) {}
