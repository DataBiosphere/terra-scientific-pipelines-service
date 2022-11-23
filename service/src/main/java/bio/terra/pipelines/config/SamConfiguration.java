package bio.terra.pipelines.config;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "pipelines.sam")
public record SamConfiguration(String basePath) {}
