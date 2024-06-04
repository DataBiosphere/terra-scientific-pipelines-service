package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

/** configuration for all properties related to cbas */
@ConfigurationProperties(prefix = "cbas")
public class CbasConfiguration {
  private Boolean callCache;
  private Boolean debugApiLogging;

  public Boolean getCallCache() {
    return callCache;
  }

  public Boolean getDebugApiLogging() {
    return debugApiLogging;
  }

  public void setCallCache(Boolean callCache) {
    this.callCache = callCache;
  }

  public void setDebugApiLogging(Boolean debugApiLogging) {
    this.debugApiLogging = debugApiLogging;
  }
}
