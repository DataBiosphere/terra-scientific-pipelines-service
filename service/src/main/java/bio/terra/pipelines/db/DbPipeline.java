package bio.terra.pipelines.db;

/** Record to hold a Pipeline record when processing in the PipelinesDao */
public record DbPipeline(String pipelineId, String displayName, String description) {}
