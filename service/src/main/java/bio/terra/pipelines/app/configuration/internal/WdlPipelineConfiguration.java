package bio.terra.pipelines.app.configuration.internal;

import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "pipelines.wdl")
@Getter
@Setter
public class WdlPipelineConfiguration {
  private Long quotaConsumedPollingIntervalSeconds;
  private boolean quotaConsumedUseCallCaching;
}
