package bio.terra.pipelines.app.configuration;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "javatemplate.sam")
public record SamConfiguration(String basePath) {}
