package bio.terra.pipelines.app.configuration.internal;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "pipelines.common")
public record PipelinesCommonConfiguration(
    Long quotaConsumedPollingIntervalSeconds,
    boolean quotaConsumedUseCallCaching,
    Long storageBucketTtlDays) {}
