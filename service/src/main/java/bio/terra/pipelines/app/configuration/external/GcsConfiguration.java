package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "gcs")
public record GcsConfiguration(long signedUrlPutDurationHours, long signedUrlGetDurationHours) {}
