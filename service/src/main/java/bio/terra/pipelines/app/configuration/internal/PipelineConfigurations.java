package bio.terra.pipelines.app.configuration.internal;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import jakarta.annotation.PostConstruct;
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
@ConfigurationProperties(prefix = "pipelines.configurations")
public class PipelineConfigurations {

  private static final Pattern PIPELINE_KEY_PATTERN = Pattern.compile("^[a-z0-9_]+_v\\d+$");

  private PipelinesCommonConfiguration common;
  private Map<String, PipelineQuotaConfig> pipelineQuotas;
  private Map<String, Map<String, PipelineConfiguration>> pipelines;

  @PostConstruct
  public void validateConfiguration() {
    validatePipelineQuotas();
    validateAllPipelineConfigurations();
  }

  private void validatePipelineQuotas() {
    if (pipelineQuotas == null || pipelineQuotas.isEmpty()) {
      throw new IllegalArgumentException("Missing pipelines.configurations.pipelineQuotas");
    }

    for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
      String pipelineName = pipelineEnum.getValue();
      PipelineQuotaConfig quota = pipelineQuotas.get(pipelineName);
      if (quota == null) {
        throw new IllegalArgumentException(
            "Missing quota definition for pipeline '%s'".formatted(pipelineName));
      }
      validateQuota(quota, pipelineName);
    }
  }

  private void validateAllPipelineConfigurations() {
    if (pipelines == null || pipelines.isEmpty()) {
      throw new IllegalArgumentException("Missing pipelines.configurations.pipelines");
    }
    for (Map.Entry<String, Map<String, PipelineConfiguration>> pipelineConfigurations :
        pipelines.entrySet()) {
      validateAllConfiguredVersionsOfPipeline(
          PipelinesEnum.enumFromStringValue(pipelineConfigurations.getKey()),
          pipelineConfigurations.getValue());
    }
  }

  /**
   * Validate each version of the pipeline supplied by versionedConfig.
   *
   * @param pipelineEnum pipeline enum being validated
   * @param versionedConfig map keyed by version string to pipeline configuration
   */
  private void validateAllConfiguredVersionsOfPipeline(
      PipelinesEnum pipelineEnum, Map<String, PipelineConfiguration> versionedConfig) {
    if (versionedConfig == null) {
      return;
    }

    for (Map.Entry<String, PipelineConfiguration> entry : versionedConfig.entrySet()) {
      int pipelineVersion;
      try {
        pipelineVersion = Integer.parseInt(entry.getKey());
      } catch (NumberFormatException e) {
        throw new IllegalArgumentException(
            "Invalid version '%s' for pipeline '%s'. Expected an integer version key"
                .formatted(entry.getKey(), pipelineEnum.getValue()),
            e);
      }
      PipelineConfiguration pipelineConfiguration = entry.getValue();
      String pipelineKey = PipelinesEnum.buildPipelineKey(pipelineEnum, pipelineVersion);

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
      String pipelineKey, PipelineConfiguration pipelineConfiguration) {
    PipelineMetadataConfig metadata = pipelineConfiguration.getMetadata();
    if (metadata == null) {
      throw new IllegalArgumentException(
          "Missing metadata for pipeline '%s'".formatted(pipelineKey));
    }

    requireText(metadata.getPipelineName(), "metadata.pipelineName", pipelineKey);
    Integer pipelineVersion =
        requireNonNull(metadata.getPipelineVersion(), "metadata.pipelineVersion", pipelineKey);
    requireText(metadata.getDisplayName(), "metadata.displayName", pipelineKey);

    String expectedKey =
        PipelinesEnum.buildPipelineKey(
            PipelinesEnum.enumFromPipelineKey(pipelineKey), pipelineVersion);
    if (!expectedKey.equals(pipelineKey)) {
      throw new IllegalArgumentException(
          "Pipeline metadata mismatch for key '%s'. Expected metadata to resolve to '%s'"
              .formatted(pipelineKey, expectedKey));
    }

    validateInputDefinitions(pipelineConfiguration.getInputDefinitionConfigs(), pipelineKey);
    validateOutputDefinitions(pipelineConfiguration.getOutputDefinitionConfigs(), pipelineKey);
  }

  private void validateInputDefinitions(
      List<PipelineInputDefinitionConfig> inputs, String pipelineKey) {
    if (inputs == null || inputs.isEmpty()) {
      throw new IllegalArgumentException(
          "Missing input definitions for pipeline '%s'".formatted(pipelineKey));
    }

    Set<String> inputNames = new HashSet<>();
    for (PipelineInputDefinitionConfig input : inputs) {
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
      List<PipelineOutputDefinitionConfig> outputs, String pipelineKey) {
    if (outputs == null || outputs.isEmpty()) {
      throw new IllegalArgumentException(
          "Missing output definitions for pipeline '%s'".formatted(pipelineKey));
    }

    Set<String> outputNames = new HashSet<>();
    for (PipelineOutputDefinitionConfig output : outputs) {
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

  private void validateQuota(PipelineQuotaConfig quota, String pipelineKey) {
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

  public PipelineQuotaConfig getQuotaForPipeline(PipelinesEnum pipelineName) {
    if (pipelineQuotas == null) {
      throw new IllegalArgumentException(
          "pipelines.configurations.pipelineQuotas is not configured");
    }
    PipelineQuotaConfig quota = pipelineQuotas.get(pipelineName.getValue());
    if (quota == null) {
      throw new IllegalArgumentException(
          "No quota configuration found for pipeline '%s'".formatted(pipelineName.getValue()));
    }
    return quota;
  }

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

  public PipelineConfiguration getPipelineConfiguration(String pipelineKey) {
    if (pipelineKey == null || !PIPELINE_KEY_PATTERN.matcher(pipelineKey).matches()) {
      throw new IllegalArgumentException(
          "Invalid pipeline key '%s'. Expected lowercase format '<pipeline_name>_v<version>'"
              .formatted(pipelineKey));
    }

    String pipelineName = PipelinesEnum.enumFromPipelineKey(pipelineKey).getValue();
    String version = String.valueOf(PipelinesEnum.versionFromPipelineKey(pipelineKey));

    Map<String, PipelineConfiguration> pipelineConfigurationMap = pipelines.get(pipelineName);
    if (pipelineConfigurationMap == null) {
      throw new IllegalArgumentException(
          "No pipeline definition map found for pipeline '%s'".formatted(pipelineName));
    }
    if (!pipelineConfigurationMap.containsKey(version)) {
      throw new IllegalArgumentException(
          "No pipeline definition found for key '%s'".formatted(pipelineKey));
    }
    return pipelines.get(pipelineName).get(version);
  }

  @Getter
  @Setter
  public static class PipelinesCommonConfiguration {
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
  public static class PipelineConfiguration {
    private PipelineMetadataConfig metadata;
    private List<PipelineInputDefinitionConfig> inputDefinitionConfigs;
    private List<PipelineOutputDefinitionConfig> outputDefinitionConfigs;

    public List<PipelineInputDefinitionConfig> getInputs() {
      return inputDefinitionConfigs;
    }

    public void setInputs(List<PipelineInputDefinitionConfig> inputs) {
      this.inputDefinitionConfigs = inputs;
    }

    public List<PipelineOutputDefinitionConfig> getOutputs() {
      return outputDefinitionConfigs;
    }

    public void setOutputs(List<PipelineOutputDefinitionConfig> outputs) {
      this.outputDefinitionConfigs = outputs;
    }
  }

  @Setter
  @Getter
  public static class PipelineMetadataConfig {
    private String pipelineName;
    private Integer pipelineVersion;
    private String displayName;
    private String description;
    private String pipelineType;
    private String toolName;
    private BigDecimal memoryRetryMultiplier;
  }

  @Setter
  @Getter
  public static class PipelineInputDefinitionConfig {
    private String name;
    private String wdlVariableName;
    private String displayName;
    private String description;
    private PipelineVariableTypesEnum type;
    private Boolean isRequired;
    private Boolean userProvided;
    private Boolean expectsCustomValue;
    private String defaultValue;
    private Double minValue;
    private Double maxValue;
    private String fileSuffix;
  }

  @Setter
  @Getter
  public static class PipelineOutputDefinitionConfig {
    private String name;
    private String wdlVariableName;
    private String displayName;
    private String description;
    private PipelineVariableTypesEnum type;
    private Boolean isRequired;
  }

  @Setter
  @Getter
  public static class PipelineQuotaConfig {
    private Integer defaultQuota;
    private Integer minQuotaConsumed;
    private QuotaUnitsEnum quotaUnits;
  }
}
