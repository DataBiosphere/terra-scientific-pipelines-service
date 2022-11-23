package bio.terra.pipelines.config;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "pipelines.status-check")
public record StatusCheckConfiguration(
    boolean enabled,
    int pollingIntervalSeconds,
    int startupWaitSeconds,
    int stalenessThresholdSeconds) {}
