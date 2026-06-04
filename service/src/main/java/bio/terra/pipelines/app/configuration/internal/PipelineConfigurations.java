package bio.terra.pipelines.app.configuration.internal;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.enumFromPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.versionFromPipelineKey;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import java.math.BigDecimal;
import java.util.List;
import java.util.Map;
import lombok.Getter;
import lombok.NoArgsConstructor;
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
@ConfigurationProperties(prefix = "configurations")
public class PipelineConfigurations {

  private CommonConfiguration common;
  private Map<String, PipelineQuotaConfiguration> pipelineQuotas;
  private Map<String, Map<String, WdlBasedPipelineConfiguration>> pipelines;

  @Getter
  @Setter
  public static class CommonConfiguration {
    private Long userDataTtlDays;
    private Long quotaConsumedPollingIntervalSeconds;
    private boolean quotaConsumedUseCallCaching;
    private Long inputQcPollingIntervalSeconds;
    private boolean inputQcUseCallCaching;
    private Long mainToolPollingIntervalSeconds;
    private boolean mainToolUseCallCaching;
    private boolean mainToolDeleteIntermediateFiles;
    private String monitoringScriptPath;
  }

  @Setter
  @Getter
  @NoArgsConstructor
  public static class WdlBasedPipelineConfiguration {
    private String displayName;
    private String description;
    private String pipelineType;
    private String toolName;
    private BigDecimal memoryRetryMultiplier;
    private List<PipelineInputDefinitionConfiguration> inputDefinitions;
    private List<PipelineOutputDefinitionConfiguration> outputDefinitions;
  }

  @Setter
  @Getter
  public static class PipelineInputDefinitionConfiguration {
    private String name;
    private String wdlVariableName;
    private String displayName;
    private String description;
    private PipelineVariableTypesEnum type;
    private Boolean isRequired;
    private Boolean userProvided;
    private String defaultValue;
    private Double minValue;
    private Double maxValue;
    private String fileSuffix;
  }

  @Setter
  @Getter
  public static class PipelineOutputDefinitionConfiguration {
    private String name;
    private String wdlVariableName;
    private String displayName;
    private String description;
    private PipelineVariableTypesEnum type;
    private Boolean isRequired;
  }

  @Setter
  @Getter
  public static class PipelineQuotaConfiguration {
    private Integer defaultQuota;
    private Integer minQuotaConsumed;
    private QuotaUnitsEnum quotaUnits;
  }

  public PipelineQuotaConfiguration getQuotaForPipeline(PipelinesEnum pipelineName) {
    if (pipelineQuotas == null) {
      throw new IllegalArgumentException(
          "pipelines.configurations.pipelineQuotas is not configured");
    }
    PipelineQuotaConfiguration quota = pipelineQuotas.get(pipelineName.getConfigKeyValue());
    if (quota == null) {
      throw new IllegalArgumentException(
          "No quota configuration found for pipeline '%s'"
              .formatted(pipelineName.getConfigKeyValue()));
    }
    return quota;
  }

  public WdlBasedPipelineConfiguration getPipelineConfiguration(String pipelineKey) {
    PipelinesEnum pipelineEnum = enumFromPipelineKey(pipelineKey);
    int pipelineVersion = versionFromPipelineKey(pipelineKey);
    String pipelineConfigKey = pipelineEnum.getConfigKeyValue();
    String version = String.valueOf(pipelineVersion);

    Map<String, WdlBasedPipelineConfiguration> pipelineConfigurationMap =
        pipelines.get(pipelineConfigKey);
    if (!pipelineConfigurationMap.containsKey(version)) {
      throw new IllegalArgumentException(
          "No pipeline definition found for key '%s'".formatted(pipelineKey));
    }

    return pipelineConfigurationMap.get(version);
  }

  public List<String> getPipelineKeys(PipelinesEnum pipelineName) {
    String pipelineConfigKey = pipelineName.getConfigKeyValue();
    Map<String, WdlBasedPipelineConfiguration> pipelineConfigurationMap =
        pipelines.get(pipelineConfigKey);

    return pipelineConfigurationMap.keySet().stream()
        .map(version -> buildPipelineKey(pipelineName, Integer.parseInt(version)))
        .sorted()
        .toList();
  }
}
