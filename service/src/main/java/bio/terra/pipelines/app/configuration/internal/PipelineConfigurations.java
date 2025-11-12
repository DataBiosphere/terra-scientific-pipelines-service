package bio.terra.pipelines.app.configuration.internal;

import java.math.BigDecimal;
import java.util.List;
import java.util.Map;
import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

/**
 * This class represents the configuration properties for pipelines. Each pipeline shares a
 * configuration class across all of its versions with each version having its own specific
 * configuration values.
 *
 * <p>This class also includes common configuration properties that apply to all pipelines.
 */
@Setter
@Getter
@ConfigurationProperties(prefix = "pipelines.configurations")
public class PipelineConfigurations {

  private PipelinesCommonConfiguration common;
  // this is a map of array imputation pipeline versions to their configurations
  private Map<String, ArrayImputationConfig> arrayImputation;

  @Getter
  @Setter
  public static class PipelinesCommonConfiguration {
    private Long userDataTtlDays;
    private Long quotaConsumedPollingIntervalSeconds;
    private boolean quotaConsumedUseCallCaching;
    private Long inputQcPollingIntervalSeconds;
    private boolean inputQcUseCallCaching;
    private String monitoringScriptPath;
  }

  @Setter
  @Getter
  public static class ArrayImputationConfig {
    private Long cromwellSubmissionPollingIntervalInSeconds;
    private boolean useCallCaching;
    private boolean deleteIntermediateFiles;
    private BigDecimal memoryRetryMultiplier;
    private Map<String, String> inputsWithCustomValues;
    private String storageWorkspaceContainerUrl;
    private List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl;

    public ArrayImputationConfig(
        Long cromwellSubmissionPollingIntervalInSeconds,
        List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl,
        String storageWorkspaceContainerUrl,
        Map<String, String> inputsWithCustomValues,
        boolean useCallCaching,
        boolean deleteIntermediateFiles,
        BigDecimal memoryRetryMultiplier) {
      this.cromwellSubmissionPollingIntervalInSeconds = cromwellSubmissionPollingIntervalInSeconds;
      this.inputKeysToPrependWithStorageWorkspaceContainerUrl =
          inputKeysToPrependWithStorageWorkspaceContainerUrl;
      this.storageWorkspaceContainerUrl = storageWorkspaceContainerUrl;

      for (Map.Entry<String, String> entry : inputsWithCustomValues.entrySet()) {
        if (entry.getValue() == null || entry.getValue().isBlank()) {
          throw new IllegalArgumentException(
              "All fields in inputsWithCustomValues must be defined. Missing value for %s"
                  .formatted(entry.getKey()));
        }
      }

      this.inputsWithCustomValues = inputsWithCustomValues;
      this.useCallCaching = useCallCaching;
      this.deleteIntermediateFiles = deleteIntermediateFiles;
      this.memoryRetryMultiplier = memoryRetryMultiplier;
    }
  }
}
