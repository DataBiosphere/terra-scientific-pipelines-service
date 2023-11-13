package bio.terra.pipelines.app.configuration;

import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Configuration;

@Configuration
@EnableConfigurationProperties
@ConfigurationProperties(prefix = "pipelines.stairway-job")
public class StairwayJobConfiguration {
  /** Number of threads to keep available */
  private int maxThreads;
  /** Timeout in seconds */
  private int timeoutSeconds;
  /** Polling interval in seconds */
  private int pollingIntervalSeconds;
  /** For identifying the application to SAM */
  private String resourceId;

  public int getTimeoutSeconds() {
    return timeoutSeconds;
  }

  public void setTimeoutSeconds(int timeoutSeconds) {
    this.timeoutSeconds = timeoutSeconds;
  }

  public int getPollingIntervalSeconds() {
    return pollingIntervalSeconds;
  }

  public void setPollingIntervalSeconds(int pollingIntervalSeconds) {
    this.pollingIntervalSeconds = pollingIntervalSeconds;
  }

  public int getMaxThreads() {
    return maxThreads;
  }

  public void setMaxThreads(int maxThreads) {
    this.maxThreads = maxThreads;
  }

  public String getResourceId() {
    return resourceId;
  }

  public void setResourceId(String resourceId) {
    this.resourceId = resourceId;
  }
}
