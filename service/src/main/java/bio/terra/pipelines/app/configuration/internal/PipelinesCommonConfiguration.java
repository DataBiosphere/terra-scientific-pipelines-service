package bio.terra.pipelines.app.configuration.internal;

import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

@Getter
@Setter
@ConfigurationProperties(prefix = "pipelines.common")
public class PipelinesCommonConfiguration {
  private Long userDataTtlDays;
  private Long quotaConsumedPollingIntervalSeconds;
  private boolean quotaConsumedUseCallCaching;
  private Long inputQcPollingIntervalSeconds;
  private boolean inputQcUseCallCaching;
}
