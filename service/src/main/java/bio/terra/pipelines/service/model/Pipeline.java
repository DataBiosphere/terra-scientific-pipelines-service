package bio.terra.pipelines.service.model;

import bio.terra.pipelines.db.DbPipeline;
import java.util.StringJoiner;

public class Pipeline {
  private final String pipelineId;
  private final String displayName;
  private final String description;

  public Pipeline(String pipelineId, String displayName, String description) {
    this.pipelineId = pipelineId;
    this.displayName = displayName;
    this.description = description;
  }

  public String getPipelineId() {
    return pipelineId;
  }

  public String getDisplayName() {
    return displayName;
  }

  public String getDescription() {
    return description;
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", Pipeline.class.getSimpleName() + "[", "]")
        .add("pipelineId=" + pipelineId)
        .add("displayName=" + displayName)
        .add("description=" + description)
        .toString();
  }

  public static Pipeline fromDb(DbPipeline dbPipeline) {
    return new Pipeline(
        dbPipeline.pipelineId(), dbPipeline.displayName(), dbPipeline.description());
  }
}
