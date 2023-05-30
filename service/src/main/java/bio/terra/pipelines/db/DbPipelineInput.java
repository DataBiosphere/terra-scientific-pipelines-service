package bio.terra.pipelines.db;

import java.util.UUID;

/** Record to hold a PipelinesInput record when processing in the JobDao */
public record DbPipelineInput(UUID jobId, String pipelineInputs) {}
