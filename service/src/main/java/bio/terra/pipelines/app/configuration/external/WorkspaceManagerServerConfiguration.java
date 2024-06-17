package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "workspace")
public record WorkspaceManagerServerConfiguration(
    String baseUri, Long sasExpirationDurationHours, Boolean debugApiLogging) {}
