package bio.terra.pipelines.app.configuration.internal;

import java.util.List;
import java.util.Map;
import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

/** configuration for properties related to imputation */
@ConfigurationProperties(prefix = "imputation")
@Getter
@Setter
public class ImputationConfiguration {

  private Long cromwellSubmissionPollingIntervalInSeconds;
  private List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl;
  private String storageWorkspaceContainerUrl;
  private Map<String, String> inputsWithCustomValues;
  private boolean useCallCaching;
  private boolean deleteIntermediateFiles;
  private boolean useReferenceDisk;

  public ImputationConfiguration(
      Long cromwellSubmissionPollingIntervalInSeconds,
      List<String> inputKeysToPrependWithStorageWorkspaceContainerUrl,
      String storageWorkspaceContainerUrl,
      Map<String, String> inputsWithCustomValues,
      boolean useCallCaching,
      boolean deleteIntermediateFiles,
      boolean useReferenceDisk) {
    this.cromwellSubmissionPollingIntervalInSeconds = cromwellSubmissionPollingIntervalInSeconds;
    this.inputKeysToPrependWithStorageWorkspaceContainerUrl =
        inputKeysToPrependWithStorageWorkspaceContainerUrl;

    if (storageWorkspaceContainerUrl.endsWith("/"))
      throw new IllegalArgumentException("Storage workspace URL must not end with a slash");

    this.storageWorkspaceContainerUrl = storageWorkspaceContainerUrl;

    for (Map.Entry<String, String> entry : inputsWithCustomValues.entrySet()) {
      if (entry.getValue() == null) {
        throw new IllegalArgumentException(
            "All fields in inputsWithCustomValues must be defined. Missing value for %s"
                .formatted(entry.getKey()));
      }
    }

    this.inputsWithCustomValues = inputsWithCustomValues;
    this.useCallCaching = useCallCaching;
    this.deleteIntermediateFiles = deleteIntermediateFiles;
    this.useReferenceDisk = useReferenceDisk;
  }
}
