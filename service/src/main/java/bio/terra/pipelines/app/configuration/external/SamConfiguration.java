package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "pipelines.sam")
public record SamConfiguration(String basePath) {}
