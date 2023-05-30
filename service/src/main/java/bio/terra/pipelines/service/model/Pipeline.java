package bio.terra.pipelines.service.model;

import bio.terra.pipelines.db.DbPipeline;
import java.util.StringJoiner;

public record Pipeline(String pipelineId, String displayName, String description) {

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
