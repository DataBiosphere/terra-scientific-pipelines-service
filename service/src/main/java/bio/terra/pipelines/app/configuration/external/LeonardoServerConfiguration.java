package bio.terra.pipelines.app.configuration.external;

import java.time.Duration;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "leonardo")
public class LeonardoServerConfiguration {

  private static final Logger log = LoggerFactory.getLogger(LeonardoServerConfiguration.class);

  private final String baseUri;
  private final List<AppType> wdsAppTypeNames;
  private final List<AppType> cbasAppTypeNames;
  private final Duration dependencyUrlCacheTtl;
  private final Boolean debugApiLogging;

  public LeonardoServerConfiguration(
      String baseUri,
      List<String> wdsAppTypeNames,
      List<String> cbasAppTypeNames,
      long dependencyUrlCacheTtlSeconds,
      Boolean debugApiLogging) {
    this.baseUri = baseUri;
    this.wdsAppTypeNames = wdsAppTypeNames.stream().map(AppType::fromValue).toList();
    this.cbasAppTypeNames = cbasAppTypeNames.stream().map(AppType::fromValue).toList();
    this.dependencyUrlCacheTtl = Duration.ofSeconds(dependencyUrlCacheTtlSeconds);
    this.debugApiLogging = debugApiLogging;
    log.info("Setting wdsAppTypes={}", wdsAppTypeNames);
    log.info("Setting cbasAppTypeNames={}", cbasAppTypeNames);
  }

  public String baseUri() {
    return baseUri;
  }

  public List<AppType> wdsAppTypeNames() {
    return wdsAppTypeNames;
  }

  public List<AppType> cbasAppTypeNames() {
    return cbasAppTypeNames;
  }

  public Duration dependencyUrlCacheTtl() {
    return dependencyUrlCacheTtl;
  }

  public Boolean debugApiLogging() {
    return debugApiLogging;
  }
}
