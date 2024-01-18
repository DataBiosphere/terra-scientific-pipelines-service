package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "cbas")
public record CbasConfiguration(Boolean debugApiLogging) {}
