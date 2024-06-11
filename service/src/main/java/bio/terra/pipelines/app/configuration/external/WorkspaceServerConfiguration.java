package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "workspace")
public record WorkspaceServerConfiguration(String baseUri, Boolean debugApiLogging) {}
