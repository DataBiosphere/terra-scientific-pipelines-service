package bio.terra.pipelines.app.configuration.internal;

import org.springframework.boot.context.properties.ConfigurationProperties;

/** configuration for all properties related to imputation */
@ConfigurationProperties(prefix = "imputation")
public class ImputationConfiguration {
  private Long cromwellSubmissionPollingIntervalInSeconds;

  public Long getCromwellSubmissionPollingIntervalInSeconds() {
    return cromwellSubmissionPollingIntervalInSeconds;
  }

  public void setCromwellSubmissionPollingIntervalInSeconds(
      Long cromwellSubmissionPollingIntervalInSeconds) {
    this.cromwellSubmissionPollingIntervalInSeconds = cromwellSubmissionPollingIntervalInSeconds;
  }
}
