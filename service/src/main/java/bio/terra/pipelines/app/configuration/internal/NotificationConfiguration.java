package bio.terra.pipelines.app.configuration.internal;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties("pipelines.notifications")
public record NotificationConfiguration(String projectId, String topicId) {}
