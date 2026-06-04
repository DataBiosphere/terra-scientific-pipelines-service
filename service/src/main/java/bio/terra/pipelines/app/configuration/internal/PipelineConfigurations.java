package bio.terra.pipelines.app.configuration.internal;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.enumFromPipelineKey;
import static bio.terra.pipelines.common.utils.PipelineKeyUtils.versionFromPipelineKey;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import java.math.BigDecimal;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.regex.Pattern;
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

  private static final Pattern PIPELINE_KEY_PATTERN = Pattern.compile("^[a-z0-9_]+_v\\d+$");

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
    if (pipelineConfigurationMap == null) {
      throw new IllegalArgumentException(
          "No pipeline definition map found for pipeline '%s'".formatted(pipelineConfigKey));
    }
    if (!pipelineConfigurationMap.containsKey(version)) {
      throw new IllegalArgumentException(
          "No pipeline definition found for key '%s'".formatted(pipelineKey));
    }

    return pipelineConfigurationMap.get(version);
  }

  // VALIDATION METHODS

  /** Comprehensive schema/content validation intended for tests that read canonical YAML. */
  public void validateConfiguration() {
    validatePipelineQuotas();
    validateAllPipelineConfigurations();
  }

  private void validatePipelineQuotas() {
    if (pipelineQuotas == null || pipelineQuotas.isEmpty()) {
      throw new IllegalArgumentException("Missing pipelines.configurations.pipelineQuotas");
    }

    for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
      String pipelineName = pipelineEnum.getConfigKeyValue();
      PipelineQuotaConfiguration quota = pipelineQuotas.get(pipelineName);
      validateQuota(quota, pipelineName);
    }
  }

  private void validateAllPipelineConfigurations() {
    if (pipelines == null || pipelines.isEmpty()) {
      throw new IllegalArgumentException("Missing pipelines.configurations.pipelines");
    }
    for (Map.Entry<String, Map<String, WdlBasedPipelineConfiguration>> pipelineConfigurations :
        pipelines.entrySet()) {
      validateAllConfiguredVersionsOfPipeline(
          PipelinesEnum.enumFromConfigKeyValue(pipelineConfigurations.getKey()),
          pipelineConfigurations.getValue());
    }
  }

  private void validateAllConfiguredVersionsOfPipeline(
      PipelinesEnum pipelineEnum, Map<String, WdlBasedPipelineConfiguration> versionedConfig) {
    if (versionedConfig == null) {
      return;
    }

    for (Map.Entry<String, WdlBasedPipelineConfiguration> entry : versionedConfig.entrySet()) {
      int pipelineVersion;
      try {
        pipelineVersion = Integer.parseInt(entry.getKey());
      } catch (NumberFormatException e) {
        throw new IllegalArgumentException(
            "Invalid version '%s' for pipeline '%s'. Expected an integer version key"
                .formatted(entry.getKey(), pipelineEnum.getConfigKeyValue()),
            e);
      }
      WdlBasedPipelineConfiguration pipelineConfiguration = entry.getValue();
      String pipelineKey = buildPipelineKey(pipelineEnum, pipelineVersion);

      if (!PIPELINE_KEY_PATTERN.matcher(pipelineKey).matches()) {
        throw new IllegalArgumentException(
            "Invalid pipeline key '%s'. Expected lowercase format '<pipeline_name>_v<version>'"
                .formatted(pipelineKey));
      }

      if (pipelineConfiguration == null) {
        throw new IllegalArgumentException(
            "Missing pipeline definition for key '%s'".formatted(pipelineKey));
      }

      validatePipelineDefinition(pipelineKey, pipelineConfiguration);
    }
  }

  private void validatePipelineDefinition(
      String pipelineKey, WdlBasedPipelineConfiguration pipelineConfiguration) {
    requireText(pipelineConfiguration.getDisplayName(), "displayName", pipelineKey);
    requireText(pipelineConfiguration.getPipelineType(), "pipelineType", pipelineKey);
    requireText(pipelineConfiguration.getToolName(), "toolName", pipelineKey);

    validateInputDefinitions(pipelineConfiguration.getInputDefinitions(), pipelineKey);
    validateOutputDefinitions(pipelineConfiguration.getOutputDefinitions(), pipelineKey);
  }

  private void validateInputDefinitions(
      List<PipelineInputDefinitionConfiguration> inputs, String pipelineKey) {
    if (inputs == null || inputs.isEmpty()) {
      throw new IllegalArgumentException(
          "Missing input definitions for pipeline '%s'".formatted(pipelineKey));
    }

    Set<String> inputNames = new HashSet<>();
    for (PipelineInputDefinitionConfiguration input : inputs) {
      if (input == null) {
        throw new IllegalArgumentException(
            "Input definition cannot be null for pipeline '%s'".formatted(pipelineKey));
      }

      String inputName = requireText(input.getName(), "inputs.name", pipelineKey);
      requireText(input.getWdlVariableName(), "inputs.wdlVariableName", pipelineKey);
      requireNonNull(input.getType(), "inputs.type", pipelineKey);
      requireNonNull(input.getIsRequired(), "inputs.isRequired", pipelineKey);
      requireNonNull(input.getUserProvided(), "inputs.userProvided", pipelineKey);

      if (!inputNames.add(inputName)) {
        throw new IllegalArgumentException(
            "Duplicate input definition name '%s' for pipeline '%s'"
                .formatted(inputName, pipelineKey));
      }
    }
  }

  private void validateOutputDefinitions(
      List<PipelineOutputDefinitionConfiguration> outputs, String pipelineKey) {
    if (outputs == null || outputs.isEmpty()) {
      throw new IllegalArgumentException(
          "Missing output definitions for pipeline '%s'".formatted(pipelineKey));
    }

    Set<String> outputNames = new HashSet<>();
    for (PipelineOutputDefinitionConfiguration output : outputs) {
      if (output == null) {
        throw new IllegalArgumentException(
            "Output definition cannot be null for pipeline '%s'".formatted(pipelineKey));
      }

      String outputName = requireText(output.getName(), "outputs.name", pipelineKey);
      requireText(output.getWdlVariableName(), "outputs.wdlVariableName", pipelineKey);
      requireNonNull(output.getType(), "outputs.type", pipelineKey);
      requireNonNull(output.getIsRequired(), "outputs.isRequired", pipelineKey);

      if (!outputNames.add(outputName)) {
        throw new IllegalArgumentException(
            "Duplicate output definition name '%s' for pipeline '%s'"
                .formatted(outputName, pipelineKey));
      }
    }
  }

  private void validateQuota(PipelineQuotaConfiguration quota, String pipelineKey) {
    if (quota == null) {
      throw new IllegalArgumentException(
          "Missing quota definition for pipeline '%s'".formatted(pipelineKey));
    }

    Integer defaultQuota =
        requireNonNull(quota.getDefaultQuota(), "quota.defaultQuota", pipelineKey);
    Integer minQuotaConsumed =
        requireNonNull(quota.getMinQuotaConsumed(), "quota.minQuotaConsumed", pipelineKey);
    requireNonNull(quota.getQuotaUnits(), "quota.quotaUnits", pipelineKey);

    if (defaultQuota < 0) {
      throw new IllegalArgumentException(
          "quota.defaultQuota must be >= 0 for pipeline '%s'".formatted(pipelineKey));
    }
    if (minQuotaConsumed < 0) {
      throw new IllegalArgumentException(
          "quota.minQuotaConsumed must be >= 0 for pipeline '%s'".formatted(pipelineKey));
    }
  }

  // TODO revisit these, consider using variable type-based validations
  private String requireText(String value, String fieldName, String pipelineKey) {
    if (value == null || value.isBlank()) {
      throw new IllegalArgumentException(
          "Missing required field '%s' for pipeline '%s'".formatted(fieldName, pipelineKey));
    }
    return value;
  }

  private <T> T requireNonNull(T value, String fieldName, String pipelineKey) {
    if (Objects.isNull(value)) {
      throw new IllegalArgumentException(
          "Missing required field '%s' for pipeline '%s'".formatted(fieldName, pipelineKey));
    }
    return value;
  }
}
