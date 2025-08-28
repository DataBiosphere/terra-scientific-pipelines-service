package bio.terra.pipelines.app.configuration.internal;

import java.math.BigDecimal;
import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

@Getter
@Setter
@ConfigurationProperties(prefix = "pipelines.common")
public class PipelinesCommonConfiguration {
  private Long quotaConsumedPollingIntervalSeconds;
  private boolean quotaConsumedUseCallCaching;
  private Long inputQcPollingIntervalSeconds;
  private boolean inputQcUseCallCaching;
  private Long userDataTtlDays;
  private BigDecimal memoryRetryMultiplier;
}
