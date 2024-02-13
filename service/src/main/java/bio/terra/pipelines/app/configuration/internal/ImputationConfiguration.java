package bio.terra.pipelines.app.configuration.internal;

import org.springframework.boot.context.properties.ConfigurationProperties;

/** configuration for all properties related to imputation */
@ConfigurationProperties(prefix = "imputation")
public class ImputationConfiguration {
  private String workspaceId;
  private Long cromwellSubmissionPollingIntervalInSeconds;

  public String getWorkspaceId() {
    return workspaceId;
  }

  public void setWorkspaceId(String workspaceId) {
    this.workspaceId = workspaceId;
  }

  public Long getCromwellSubmissionPollingIntervalInSeconds() {
    return cromwellSubmissionPollingIntervalInSeconds;
  }

  public void setCromwellSubmissionPollingIntervalInSeconds(
      Long cromwellSubmissionPollingIntervalInSeconds) {
    this.cromwellSubmissionPollingIntervalInSeconds = cromwellSubmissionPollingIntervalInSeconds;
  }
}
