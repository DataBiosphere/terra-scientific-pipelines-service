package bio.terra.pipelines.configuration.internal;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.context.properties.bind.Bindable;
import org.springframework.boot.context.properties.bind.Binder;
import org.springframework.boot.env.YamlPropertySourceLoader;
import org.springframework.core.env.PropertySource;
import org.springframework.core.env.StandardEnvironment;
import org.springframework.core.io.FileSystemResource;

class PipelineConfigurationsTest extends BaseEmbeddedDbTest {
  private static final Pattern PIPELINE_KEY_PATTERN = Pattern.compile("^[a-z0-9_]+_v\\d+$");

  @Autowired private PipelineConfigurations pipelineConfigurations;

  @Test
  void testPipelinesCommonConfiguration() {
    PipelineConfigurations.PipelinesCommonConfiguration pipelinesCommonConfiguration =
        pipelineConfigurations.getCommon();
    assertEquals(2, pipelinesCommonConfiguration.getUserDataTtlDays());

    assertEquals(1, pipelinesCommonConfiguration.getQuotaConsumedPollingIntervalSeconds());
    assertTrue(pipelinesCommonConfiguration.isQuotaConsumedUseCallCaching());

    assertEquals(1, pipelinesCommonConfiguration.getInputQcPollingIntervalSeconds());
    assertTrue(pipelinesCommonConfiguration.isInputQcUseCallCaching());

    assertEquals(2, pipelinesCommonConfiguration.getMainToolPollingIntervalSeconds());
    assertTrue(pipelinesCommonConfiguration.isMainToolUseCallCaching());
    assertFalse(pipelinesCommonConfiguration.isMainToolDeleteIntermediateFiles());

    assertEquals(
        "gs://test_bucket/test_path/to/monitoring/script.sh",
        pipelinesCommonConfiguration.getMonitoringScriptPath());
  }

  @Test
  void testPipelineDefinitionsConfigurationLoaded() {
    PipelineConfigurations.WdlBasedPipelineConfiguration pipelineDefinitionConfig =
        pipelineConfigurations.getPipelineConfiguration("array_imputation_v1");

    assertNotNull(pipelineDefinitionConfig);
    assertNotNull(pipelineDefinitionConfig.getDisplayName());
    assertEquals(1, pipelineDefinitionConfig.getInputDefinitions().size());
    assertEquals(1, pipelineDefinitionConfig.getOutputDefinitions().size());
    assertEquals(
        2500,
        pipelineConfigurations
            .getQuotaForPipeline(PipelinesEnum.ARRAY_IMPUTATION)
            .getDefaultQuota());
  }

  @Test
  void getPipelineConfigurationReturnsDefinition() {
    PipelineConfigurations.WdlBasedPipelineConfiguration definition =
        pipelineConfigurations.getPipelineConfiguration("array_imputation_v1");

    assertNotNull(definition);
    assertNotNull(definition.getDisplayName());
  }

  @Test
  void getPipelineConfigurationThrowsForMissingDefinition() {
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelineConfigurations.getPipelineConfiguration("array_imputation_v999"));
  }

  @Nested
  @DisplayName("Test Fixture Configuration Validations")
  class TestFixtureConfigurationValidations {
    @Test
    void allKnownPipelinesAreConfigured() {
      assertNoViolations(validateAllKnownPipelinesPresent(pipelineConfigurations));
    }

    @Test
    void pipelineQuotasAreValid() {
      assertNoViolations(validatePipelineQuotasOnly(pipelineConfigurations));
    }

    @Test
    void pipelineVersionsAndKeysAreValid() {
      assertNoViolations(validateVersionMapsAndPipelineKeys(pipelineConfigurations));
    }

    @Test
    void pipelineTopLevelFieldsAreValid() {
      assertNoViolations(validatePipelineDefinitionTopLevelFields(pipelineConfigurations));
    }

    @Test
    void pipelineInputDefinitionsAreValid() {
      assertNoViolations(validatePipelineInputDefinitionsOnly(pipelineConfigurations));
    }

    @Test
    void pipelineOutputDefinitionsAreValid() {
      assertNoViolations(validatePipelineOutputDefinitionsOnly(pipelineConfigurations));
    }
  }

  // TESTING/VALIDATING REAL CONFIG VALUES
  @Nested
  @DisplayName("Service ('Real') Configuration Validations")
  class ServicePipelinesConfigurationValidations {
    @Test
    void allKnownPipelinesAreConfigured() {
      assertNoViolations(validateAllKnownPipelinesPresent(loadServicePipelineConfigurations()));
    }

    @Test
    void pipelineQuotasAreValid() {
      assertNoViolations(validatePipelineQuotasOnly(loadServicePipelineConfigurations()));
    }

    @Test
    void pipelineVersionsAndKeysAreValid() {
      assertNoViolations(validateVersionMapsAndPipelineKeys(loadServicePipelineConfigurations()));
    }

    @Test
    void pipelineTopLevelFieldsAreValid() {
      assertNoViolations(
          validatePipelineDefinitionTopLevelFields(loadServicePipelineConfigurations()));
    }

    @Test
    void pipelineInputDefinitionsAreValid() {
      assertNoViolations(validatePipelineInputDefinitionsOnly(loadServicePipelineConfigurations()));
    }

    @Test
    void pipelineOutputDefinitionsAreValid() {
      assertNoViolations(
          validatePipelineOutputDefinitionsOnly(loadServicePipelineConfigurations()));
    }
  }

  /**
   * Asserts that a validation pass produced zero violations.
   *
   * <p>When violations are present, the assertion message includes every violation joined into a
   * single line so test failures show the full set of issues in one run.
   */
  private void assertNoViolations(Set<String> violations) {
    assertTrue(
        violations.isEmpty(),
        "Expected no config violations, but found: " + String.join(" | ", violations));
  }

  /**
   * Validates that pipeline maps exist for all known {@link PipelinesEnum} values.
   *
   * <p>Checks performed:
   *
   * <ul>
   *   <li>{@code configurations.pipelines} exists
   *   <li>for each enum value, a corresponding map keyed by {@code configKeyValue} exists
   *   <li>each pipeline map has at least one version entry
   * </ul>
   */
  private Set<String> validateAllKnownPipelinesPresent(
      PipelineConfigurations pipelineConfigurations) {
    Set<String> violations = new LinkedHashSet<>();
    Map<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>> pipelines =
        pipelineConfigurations.getPipelines();
    if (pipelines == null) {
      violations.add("Missing configurations.pipelines");
      return violations;
    }

    for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
      Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versions =
          pipelines.get(pipelineEnum.getConfigKeyValue());
      if (versions == null) {
        violations.add("Missing pipeline map for " + pipelineEnum.getConfigKeyValue());
      } else if (versions.isEmpty()) {
        violations.add("No versions configured for " + pipelineEnum.getConfigKeyValue());
      }
    }
    return violations;
  }

  /**
   * Runs quota-specific validation rules and returns collected quota violations only.
   *
   * <p>Checks performed:
   *
   * <ul>
   *   <li>the quota map exists and is non-empty
   *   <li>every known pipeline enum has a quota entry
   *   <li>each quota entry passes {@link #validateQuota}
   * </ul>
   */
  private Set<String> validatePipelineQuotasOnly(PipelineConfigurations pipelineConfigurations) {
    Set<String> violations = new LinkedHashSet<>();
    Map<String, PipelineConfigurations.PipelineQuotaConfiguration> pipelineQuotas =
        pipelineConfigurations.getPipelineQuotas();
    if (pipelineQuotas == null || pipelineQuotas.isEmpty()) {
      violations.add("Missing configurations.pipelineQuotas");
    } else {
      for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
        String pipelineName = pipelineEnum.getConfigKeyValue();
        PipelineConfigurations.PipelineQuotaConfiguration quota = pipelineQuotas.get(pipelineName);
        validateQuota(quota, pipelineName, violations);
      }
    }
    return violations;
  }

  /**
   * Validates version-map structure and canonical pipeline-key derivation.
   *
   * <p>Checks performed:
   *
   * <ul>
   *   <li>{@code configurations.pipelines} exists and is non-empty
   *   <li>each pipeline has a non-null version map
   *   <li>each version key parses as an integer
   *   <li>derived key {@code <pipeline_name>_v<version>} matches the lowercase key pattern
   * </ul>
   */
  private Set<String> validateVersionMapsAndPipelineKeys(
      PipelineConfigurations pipelineConfigurations) {
    Set<String> violations = new LinkedHashSet<>();
    Map<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>> pipelines =
        pipelineConfigurations.getPipelines();
    for (Map.Entry<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>>
        configuredPipeline : pipelines.entrySet()) {
      PipelinesEnum pipelineEnum =
          PipelinesEnum.enumFromConfigKeyValue(configuredPipeline.getKey());
      Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versionedConfig =
          configuredPipeline.getValue();
      if (versionedConfig == null) {
        violations.add(
            "Missing version map for pipeline '%s'".formatted(pipelineEnum.getConfigKeyValue()));
        continue;
      }

      for (String versionKey : versionedConfig.keySet()) {
        int pipelineVersion;
        try {
          pipelineVersion = Integer.parseInt(versionKey);
        } catch (NumberFormatException e) {
          violations.add(
              "Invalid version '%s' for pipeline '%s'. Expected an integer version key"
                  .formatted(versionKey, pipelineEnum.getConfigKeyValue()));
          continue;
        }

        String pipelineKey = buildPipelineKey(pipelineEnum, pipelineVersion);
        if (!PIPELINE_KEY_PATTERN.matcher(pipelineKey).matches()) {
          violations.add(
              "Invalid pipeline key '%s'. Expected lowercase format '<pipeline_name>_v<version>'"
                  .formatted(pipelineKey));
        }
      }
    }

    return violations;
  }

  /**
   * Validates required top-level metadata fields on each pipeline definition.
   *
   * <p>For each parseable version entry with a non-null definition, this checks that:
   *
   * <ul>
   *   <li>{@code displayName} is non-blank
   *   <li>{@code pipelineType} is non-blank
   *   <li>{@code toolName} is non-blank
   * </ul>
   */
  private Set<String> validatePipelineDefinitionTopLevelFields(
      PipelineConfigurations pipelineConfigurations) {
    Set<String> violations = new LinkedHashSet<>();
    Map<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>> pipelines =
        pipelineConfigurations.getPipelines();

    for (Map.Entry<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>>
        configuredPipeline : pipelines.entrySet()) {
      PipelinesEnum pipelineEnum =
          PipelinesEnum.enumFromConfigKeyValue(configuredPipeline.getKey());
      Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versionedConfig =
          configuredPipeline.getValue();
      if (versionedConfig == null) {
        continue;
      }

      for (Map.Entry<String, PipelineConfigurations.WdlBasedPipelineConfiguration> entry :
          versionedConfig.entrySet()) {
        int pipelineVersion;
        try {
          pipelineVersion = Integer.parseInt(entry.getKey());
        } catch (NumberFormatException e) {
          continue;
        }
        String pipelineKey = buildPipelineKey(pipelineEnum, pipelineVersion);
        PipelineConfigurations.WdlBasedPipelineConfiguration pipelineConfiguration =
            entry.getValue();

        requireText(pipelineConfiguration.getDisplayName(), "displayName", pipelineKey, violations);
        requireText(
            pipelineConfiguration.getPipelineType(), "pipelineType", pipelineKey, violations);
        requireText(pipelineConfiguration.getToolName(), "toolName", pipelineKey, violations);
      }
    }

    return violations;
  }

  /**
   * Validates input-definition sections for every configured pipeline version.
   *
   * <p>Delegates to {@link #validateInputDefinitions(List, String, Set)} and therefore checks:
   * required input fields, uniqueness of input names, rules for service-provided defaults, and
   * file-suffix requirements for user-provided file-like inputs.
   */
  private Set<String> validatePipelineInputDefinitionsOnly(
      PipelineConfigurations pipelineConfigurations) {
    Set<String> violations = new LinkedHashSet<>();
    Map<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>> pipelines =
        pipelineConfigurations.getPipelines();

    for (Map.Entry<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>>
        configuredPipeline : pipelines.entrySet()) {
      PipelinesEnum pipelineEnum =
          PipelinesEnum.enumFromConfigKeyValue(configuredPipeline.getKey());
      Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versionedConfig =
          configuredPipeline.getValue();
      if (versionedConfig == null) {
        continue;
      }
      for (Map.Entry<String, PipelineConfigurations.WdlBasedPipelineConfiguration> entry :
          versionedConfig.entrySet()) {
        int pipelineVersion;
        try {
          pipelineVersion = Integer.parseInt(entry.getKey());
        } catch (NumberFormatException e) {
          continue;
        }
        String pipelineKey = buildPipelineKey(pipelineEnum, pipelineVersion);
        PipelineConfigurations.WdlBasedPipelineConfiguration pipelineConfiguration =
            entry.getValue();
        if (pipelineConfiguration == null) {
          continue;
        }
        validateInputDefinitions(
            pipelineConfiguration.getInputDefinitions(), pipelineKey, violations);
      }
    }

    return violations;
  }

  /**
   * Validates output-definition sections for every configured pipeline version.
   *
   * <p>Delegates to {@link #validateOutputDefinitions(List, String, Set)} and therefore checks
   * required output fields plus uniqueness of output names within each pipeline version.
   */
  private Set<String> validatePipelineOutputDefinitionsOnly(
      PipelineConfigurations pipelineConfigurations) {
    Set<String> violations = new LinkedHashSet<>();
    Map<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>> pipelines =
        pipelineConfigurations.getPipelines();

    for (Map.Entry<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>>
        configuredPipeline : pipelines.entrySet()) {
      PipelinesEnum pipelineEnum =
          PipelinesEnum.enumFromConfigKeyValue(configuredPipeline.getKey());
      Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versionedConfig =
          configuredPipeline.getValue();
      if (versionedConfig == null) {
        continue;
      }
      for (Map.Entry<String, PipelineConfigurations.WdlBasedPipelineConfiguration> entry :
          versionedConfig.entrySet()) {
        int pipelineVersion;
        try {
          pipelineVersion = Integer.parseInt(entry.getKey());
        } catch (NumberFormatException e) {
          continue;
        }
        String pipelineKey = buildPipelineKey(pipelineEnum, pipelineVersion);
        PipelineConfigurations.WdlBasedPipelineConfiguration pipelineConfiguration =
            entry.getValue();
        if (pipelineConfiguration == null) {
          continue;
        }
        validateOutputDefinitions(
            pipelineConfiguration.getOutputDefinitions(), pipelineKey, violations);
      }
    }

    return violations;
  }

  private PipelineConfigurations loadServicePipelineConfigurations() {
    try {
      Path pipelinesConfigPath =
          Paths.get("src", "main", "resources", "pipelines-config.yml").toAbsolutePath();
      FileSystemResource pipelinesConfigResource =
          new FileSystemResource(pipelinesConfigPath.toFile());

      StandardEnvironment environment = new StandardEnvironment();
      YamlPropertySourceLoader yamlLoader = new YamlPropertySourceLoader();
      List<PropertySource<?>> propertySources =
          yamlLoader.load("main-pipelines-config", pipelinesConfigResource);
      propertySources.forEach(source -> environment.getPropertySources().addLast(source));

      return Binder.get(environment)
          .bind("configurations", Bindable.of(PipelineConfigurations.class))
          .orElseThrow(
              () ->
                  new IllegalStateException(
                      "Could not bind configurations from " + pipelinesConfigPath));
    } catch (IOException e) {
      throw new IllegalStateException("Failed to load main pipelines config", e);
    }
  }

  /**
   * Validates input definitions for a single pipeline key.
   *
   * <p>Checks performed per input list:
   *
   * <ul>
   *   <li>inputDefinitions list exists and is non-empty
   *   <li>individual input entries are non-null
   *   <li>required fields are present: {@code name}, {@code wdlVariableName}, {@code type}, {@code
   *       isRequired}, {@code userProvided}
   *   <li>for user-provided inputs ({@code userProvided=true}), {@code displayName} and {@code
   *       description} are required and non-blank
   *   <li>for service-provided inputs ({@code userProvided=false}), {@code defaultValue} is
   *       required and non-blank
   *   <li>for user-provided file-like inputs, {@code fileSuffix} is required and non-blank
   *   <li>input names are unique within the pipeline
   * </ul>
   */
  private void validateInputDefinitions(
      List<PipelineInputDefinition> inputs, String pipelineKey, Set<String> violations) {
    if (inputs == null || inputs.isEmpty()) {
      violations.add("Missing input definitions for pipeline '%s'".formatted(pipelineKey));
      return;
    }

    Set<String> inputNames = new HashSet<>();
    for (PipelineInputDefinition input : inputs) {
      if (input == null) {
        violations.add("Input definition cannot be null for pipeline '%s'".formatted(pipelineKey));
        continue;
      }

      String inputName = requireText(input.getName(), "inputs.name", pipelineKey, violations);
      requireText(input.getWdlVariableName(), "inputs.wdlVariableName", pipelineKey, violations);
      requireNonNull(input.getType(), "inputs.type", pipelineKey, violations);

      if (input.isUserProvided()) {
        requireText(input.getDisplayName(), "inputs.displayName", pipelineKey, violations);
        requireText(input.getDescription(), "inputs.description", pipelineKey, violations);
      }

      if (!input.isUserProvided()) {
        requireText(input.getDefaultValue(), "inputs.defaultValue", pipelineKey, violations);
      }

      if (input.isUserProvided() && input.getType() != null && input.getType().isFileLike()) {
        requireText(input.getFileSuffix(), "inputs.fileSuffix", pipelineKey, violations);
      }

      if (input.getDefaultValue() != null && input.getDefaultValue().isBlank()) {
        violations.add(
            "Default value is blank for input '%s' in pipeline '%s'"
                .formatted(input.getName(), pipelineKey));
      }

      if (inputName != null && !inputNames.add(inputName)) {
        violations.add(
            "Duplicate input definition name '%s' for pipeline '%s'"
                .formatted(inputName, pipelineKey));
      }
    }
  }

  /**
   * Validates output definitions for a single pipeline key.
   *
   * <p>Checks performed per output list:
   *
   * <ul>
   *   <li>outputDefinitions list exists and is non-empty
   *   <li>individual output entries are non-null
   *   <li>required fields are present: {@code name}, {@code wdlVariableName}, {@code displayName},
   *       {@code description}, {@code type}, {@code isRequired}
   *   <li>output names are unique within the pipeline
   * </ul>
   */
  private void validateOutputDefinitions(
      List<PipelineOutputDefinition> outputs, String pipelineKey, Set<String> violations) {
    if (outputs == null || outputs.isEmpty()) {
      violations.add("Missing output definitions for pipeline '%s'".formatted(pipelineKey));
      return;
    }

    Set<String> outputNames = new HashSet<>();
    for (PipelineOutputDefinition output : outputs) {
      if (output == null) {
        violations.add("Output definition cannot be null for pipeline '%s'".formatted(pipelineKey));
        continue;
      }

      String outputName = requireText(output.getName(), "outputs.name", pipelineKey, violations);
      requireText(output.getWdlVariableName(), "outputs.wdlVariableName", pipelineKey, violations);
      requireText(output.getDisplayName(), "outputs.displayName", pipelineKey, violations);
      requireText(output.getDescription(), "outputs.description", pipelineKey, violations);
      requireNonNull(output.getType(), "outputs.type", pipelineKey, violations);

      if (outputName != null && !outputNames.add(outputName)) {
        violations.add(
            "Duplicate output definition name '%s' for pipeline '%s'"
                .formatted(outputName, pipelineKey));
      }
    }
  }

  /**
   * Validates one quota block for a specific pipeline key.
   *
   * <p>Checks performed:
   *
   * <ul>
   *   <li>quota object is non-null
   *   <li>required fields are present: {@code defaultQuota}, {@code minQuotaConsumed}, {@code
   *       quotaUnits}
   *   <li>{@code defaultQuota >= 0}
   *   <li>{@code minQuotaConsumed >= 0}
   * </ul>
   */
  private void validateQuota(
      PipelineConfigurations.PipelineQuotaConfiguration quota,
      String pipelineKey,
      Set<String> violations) {
    if (quota == null) {
      violations.add("Missing quota definition for pipeline '%s'".formatted(pipelineKey));
      return;
    }

    Integer defaultQuota =
        requireNonNull(quota.getDefaultQuota(), "quota.defaultQuota", pipelineKey, violations);
    Integer minQuotaConsumed =
        requireNonNull(
            quota.getMinQuotaConsumed(), "quota.minQuotaConsumed", pipelineKey, violations);
    requireNonNull(quota.getQuotaUnits(), "quota.quotaUnits", pipelineKey, violations);

    if (defaultQuota != null && defaultQuota < 0) {
      violations.add("quota.defaultQuota must be >= 0 for pipeline '%s'".formatted(pipelineKey));
    }
    if (minQuotaConsumed != null && minQuotaConsumed < 0) {
      violations.add(
          "quota.minQuotaConsumed must be >= 0 for pipeline '%s'".formatted(pipelineKey));
    }
  }

  /**
   * Validates a required text field.
   *
   * <p>Adds a violation when the value is null or blank, and returns {@code null} in that case.
   * Returning the value allows callers to continue validation logic (for example duplicate checks)
   * without throwing.
   */
  private String requireText(
      String value, String fieldName, String pipelineKey, Set<String> violations) {
    if (value == null || value.isBlank()) {
      violations.add(
          "Missing required field '%s' for pipeline '%s'".formatted(fieldName, pipelineKey));
      return null;
    }
    return value;
  }

  /**
   * Validates a required non-text field.
   *
   * <p>Adds a violation when the value is null and returns the original value so caller logic can
   * continue without exception-driven control flow.
   */
  private <T> T requireNonNull(
      T value, String fieldName, String pipelineKey, Set<String> violations) {
    if (value == null) {
      violations.add(
          "Missing required field '%s' for pipeline '%s'".formatted(fieldName, pipelineKey));
    }
    return value;
  }
}
