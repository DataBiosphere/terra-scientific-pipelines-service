package bio.terra.pipelines.db;

import java.sql.Timestamp;
import java.util.UUID;

/** Record to hold a Job record when processing in the JobDao */
public record DbJob(
    UUID jobId,
    String userId,
    String pipelineId,
    String pipelineVersion,
    Timestamp timeSubmitted,
    Timestamp timeCompleted,
    String status) {}
