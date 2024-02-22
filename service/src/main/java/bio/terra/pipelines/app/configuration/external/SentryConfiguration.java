package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "pipelines.sentry")
public record SentryConfiguration(String dsn, String environment) {}
