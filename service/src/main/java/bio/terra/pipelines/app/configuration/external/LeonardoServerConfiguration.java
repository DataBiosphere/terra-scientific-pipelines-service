package bio.terra.pipelines.app.configuration.external;

import java.time.Duration;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.context.properties.ConstructorBinding;

/**
 * @param baseUri - Leonardo URI to send requests to
 * @param wdsAppTypeNames - names used to signify what is the WDS app we want to use
 * @param cbasAppTypeNames - names used to signify what is the cbas app we want to use
 * @param dependencyUrlCacheTtl - how long (in seconds) to keep items in the cache
 * @param debugApiLogging
 */
@ConfigurationProperties(prefix = "leonardo")
public record LeonardoServerConfiguration(
    String baseUri,
    List<AppType> wdsAppTypeNames,
    List<AppType> cbasAppTypeNames,
    Duration dependencyUrlCacheTtl,
    Boolean debugApiLogging) {

  private static final Logger log = LoggerFactory.getLogger(LeonardoServerConfiguration.class);

  @ConstructorBinding
  public LeonardoServerConfiguration(
      String baseUri,
      List<String> wdsAppTypeNames,
      List<String> cbasAppTypeNames,
      long dependencyUrlCacheTtlSeconds,
      Boolean debugApiLogging) {
    this(
        baseUri,
        wdsAppTypeNames.stream().map(AppType::fromValue).toList(),
        cbasAppTypeNames.stream().map(AppType::fromValue).toList(),
        Duration.ofSeconds(dependencyUrlCacheTtlSeconds),
        debugApiLogging);
    log.info("Setting wdsAppTypes={}", wdsAppTypeNames);
    log.info("Setting cbasAppTypeNames={}", cbasAppTypeNames);
  }
}
