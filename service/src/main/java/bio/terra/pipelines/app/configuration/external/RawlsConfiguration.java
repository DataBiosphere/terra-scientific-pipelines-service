package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "rawls")
public record RawlsConfiguration(String baseUri, Boolean debugApiLogging) {}
