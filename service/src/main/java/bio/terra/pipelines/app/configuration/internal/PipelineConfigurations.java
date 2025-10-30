package bio.terra.pipelines.app.configuration.internal;

import java.util.Map;
import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

@Setter
@Getter
@ConfigurationProperties(prefix = "pipelines.configurations.pipelines")
public class PipelineConfigurations {

  private Map<String, ImputationConfig> array_imputation;
  private Map<String, OtherPipelineConfig> other_pipeline;

  @Setter
  @Getter
  public static class ImputationConfig {
    private int cromwellSubmissionPollingIntervalInSeconds;
    private boolean useCallCaching;
    private boolean deleteIntermediateFiles;
    private double memoryRetryMultiplier;
    private Map<String, String> inputsWithCustomValues;
    private String storageWorkspaceContainerUrl;
    private String[] inputKeysToPrependWithStorageWorkspaceContainerUrl;
  }

  @Setter
  @Getter
  public static class OtherPipelineConfig {
    private String differentField;
    private boolean otherFlag;
    // ... other pipeline-specific fields
  }
}
