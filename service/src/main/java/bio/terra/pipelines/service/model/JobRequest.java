package bio.terra.pipelines.service.model;

import com.fasterxml.jackson.databind.annotation.JsonDeserialize;

@JsonDeserialize(builder = JobRequest.Builder.class)
public class JobRequest {
  private final String pipelineId;
  private final String pipelineVersion;

  public JobRequest(String pipelineId, String pipelineVersion) {

    this.pipelineId = pipelineId;
    this.pipelineVersion = pipelineVersion;
  }

  public String getPipelineId() {
    return pipelineId;
  }

  public String getPipelineVersion() {
    return pipelineVersion;
  }

  public static class Builder {
    private String pipelineId;
    private String pipelineVersion;

    public Builder setPipelineId(String pipelineId) {
      this.pipelineId = pipelineId;
      return this;
    }

    public Builder setPipelineVersion(String pipelineVersion) {
      this.pipelineVersion = pipelineVersion;
      return this;
    }

    public JobRequest build() {
      return new JobRequest(pipelineId, pipelineVersion);
    }
  }
}
