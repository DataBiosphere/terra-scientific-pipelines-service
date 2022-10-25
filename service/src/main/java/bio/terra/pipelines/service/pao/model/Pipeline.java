package bio.terra.pipelines.service.pao.model;

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

  public String toShortString() {
    return String.format("%s:%s (%s)", pipelineId, displayName, description);
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
    return new Pipeline.Builder()
        .setPipelineId(dbPipeline.pipelineId())
        .setDisplayName(dbPipeline.displayName())
        .setDescription(dbPipeline.description())
        .build();
  }

  public static class Builder {
    private String pipelineId;
    private String displayName;
    private String description;

    public Builder setPipelineId(String pipelineId) {
      this.pipelineId = pipelineId;
      return this;
    }

    public Builder setDisplayName(String displayName) {
      this.displayName = displayName;
      return this;
    }

    public Builder setDescription(String description) {
      this.description = description;
      return this;
    }

    public Pipeline build() {
      return new Pipeline(pipelineId, displayName, description);
    }
  }
}
