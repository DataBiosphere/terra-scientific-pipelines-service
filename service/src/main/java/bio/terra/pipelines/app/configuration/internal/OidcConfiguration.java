package bio.terra.pipelines.app.configuration.internal;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties("oidc")
public record OidcConfiguration(String clientId, String authorityEndpoint) {}
