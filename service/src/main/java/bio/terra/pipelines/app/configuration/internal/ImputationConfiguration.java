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
  private List<String> inputKeysToPrependWithStorageUrl;
  private String storageWorkspaceStorageUrl;
  private Map<String, Object> inputsWithCustomValues;
  private boolean useCallCaching;
  private boolean deleteIntermediateFiles;
  private boolean useReferenceDisk;
}
