package bio.terra.pipelines.app.configuration.internal;

import java.util.List;
import org.springframework.boot.context.properties.ConfigurationProperties;

/** configuration for properties related to imputation */
@ConfigurationProperties(prefix = "imputation")
public class ImputationConfiguration {
  private Long cromwellSubmissionPollingIntervalInSeconds;
  private List<String> inputKeysToPrependWithStorageUrl;

  private String storageWorkspaceStorageUrl;

  public Long getCromwellSubmissionPollingIntervalInSeconds() {
    return cromwellSubmissionPollingIntervalInSeconds;
  }

  public List<String> getInputKeysToPrependWithStorageUrl() {
    return inputKeysToPrependWithStorageUrl;
  }

  public String getStorageWorkspaceStorageUrl() {
    return storageWorkspaceStorageUrl;
  }

  public void setCromwellSubmissionPollingIntervalInSeconds(
      Long cromwellSubmissionPollingIntervalInSeconds) {
    this.cromwellSubmissionPollingIntervalInSeconds = cromwellSubmissionPollingIntervalInSeconds;
  }

  public void setInputKeysToPrependWithStorageUrl(List<String> inputKeysToPrependWithStorageUrl) {
    this.inputKeysToPrependWithStorageUrl = inputKeysToPrependWithStorageUrl;
  }

  public void setStorageWorkspaceStorageUrl(String storageWorkspaceStorageUrl) {
    this.storageWorkspaceStorageUrl = storageWorkspaceStorageUrl;
  }
}
