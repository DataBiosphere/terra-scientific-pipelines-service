package bio.terra.pipelines.app.configuration.internal;

import java.math.BigDecimal;
import java.util.List;
import java.util.Map;
import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

@Setter
@Getter
@ConfigurationProperties(prefix = "pipelines.configurations")
public class PipelineConfigurations {

  private Map<String, ImputationConfig> arrayImputation;

  @Setter
  @Getter
  public static class ImputationConfig {
    private Long cromwellSubmissionPollingIntervalInSeconds;
    private boolean useCallCaching;
    private boolean deleteIntermediateFiles;
    private BigDecimal memoryRetryMultiplier;
    private Map<String, String> inputsWithCustomValues;
    private String storageWorkspaceContainerUrl;
    private List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl;

    public ImputationConfig(
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
