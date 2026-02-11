package bio.terra.pipelines.app.configuration.internal;

import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.context.annotation.Configuration;

@Getter
@Setter
@Configuration
@ConfigurationProperties(prefix = "pipelines.cloud-integration")
public class CloudIntegrationConfiguration {
  private String serviceAccountGroup;
}
