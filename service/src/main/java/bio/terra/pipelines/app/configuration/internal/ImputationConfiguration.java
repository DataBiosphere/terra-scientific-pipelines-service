package bio.terra.pipelines.app.configuration.internal;

import org.springframework.boot.context.properties.ConfigurationProperties;

/** configuration for all properties related to imputation */
@ConfigurationProperties(prefix = "imputation")
public record ImputationConfiguration(String workspaceId) {}
