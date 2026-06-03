package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.fasterxml.jackson.core.type.TypeReference;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;
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

  @Autowired private PipelineConfigurations pipelineConfigurations;

  @Test
  void testArrayImputationV1Configuration() {
    // note these are the values in test/resources/pipelines-config.yml and not production values
    PipelineConfigurations.PipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getPipelineConfiguration("array_imputation_v1");
    PipelineConfigurations.PipelineMetadataConfig metadata = pipelineConfiguration.getMetadata();

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(0.0);

    assertEquals(memoryRetryMultiplier, metadata.getMemoryRetryMultiplier());
  }

  @Test
  void testArrayImputationV2Configuration() {
    // note these are the values in test/resources/pipelines-config.yml and not production values
    PipelineConfigurations.PipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getPipelineConfiguration("array_imputation_v2");
    PipelineConfigurations.PipelineMetadataConfig metadata = pipelineConfiguration.getMetadata();

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(1.4);

    assertEquals(memoryRetryMultiplier, metadata.getMemoryRetryMultiplier());
  }

  @Test
  void testLowPassImputationV1Configuration() {
    // note these are the values in test/resources/pipelines-config.yml and not production values
    PipelineConfigurations.PipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getPipelineConfiguration("low_pass_imputation_v1");
    PipelineConfigurations.PipelineMetadataConfig metadata = pipelineConfiguration.getMetadata();

    assertEquals(BigDecimal.valueOf(2.0), metadata.getMemoryRetryMultiplier());
  }

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
    PipelineConfigurations.PipelineConfiguration pipelineDefinitionConfig =
        pipelineConfigurations.getPipelineConfiguration("array_imputation_v1");

    assertNotNull(pipelineDefinitionConfig);
    assertEquals("array_imputation", pipelineDefinitionConfig.getMetadata().getPipelineName());
    assertEquals(1, pipelineDefinitionConfig.getMetadata().getPipelineVersion());
    assertEquals(1, pipelineDefinitionConfig.getInputDefinitionConfigs().size());
    assertEquals(1, pipelineDefinitionConfig.getOutputDefinitionConfigs().size());
    assertEquals(
        100,
        pipelineConfigurations
            .getQuotaForPipeline(PipelinesEnum.ARRAY_IMPUTATION)
            .getDefaultQuota());
  }

  @Test
  void invalidPipelineKeyInPipelineDefinitionsThrows() {
    PipelineConfigurations pipelineConfigurationsUnderTest = new PipelineConfigurations();
    pipelineConfigurationsUnderTest.setPipelines(
        Map.of(
            "array_imputation",
            Map.of("BadVersion", buildValidTestPipelineDefinition("array_imputation", 1))));
    pipelineConfigurationsUnderTest.setPipelineQuotas(
        Map.of(
            "array_imputation",
            buildValidTestQuotaConfig(),
            "low_pass_imputation",
            buildValidTestQuotaConfig()));

    assertThrows(
        IllegalArgumentException.class, pipelineConfigurationsUnderTest::validateConfiguration);
  }

  @Test
  void missingPipelineQuotaThrows() {
    PipelineConfigurations pipelineConfigurationsUnderTest = new PipelineConfigurations();
    pipelineConfigurationsUnderTest.setPipelines(
        Map.of(
            "array_imputation",
            Map.of("1", buildValidTestPipelineDefinition("array_imputation", 1)),
            "low_pass_imputation",
            Map.of("1", buildValidTestPipelineDefinition("low_pass_imputation", 1))));
    pipelineConfigurationsUnderTest.setPipelineQuotas(
        Map.of("array_imputation", buildValidTestQuotaConfig()));

    assertThrows(
        IllegalArgumentException.class, pipelineConfigurationsUnderTest::validateConfiguration);
  }

  @Test
  void runtimeValidationMissingServiceProvidedDefaultThrows() {
    PipelineConfigurations pipelineConfigurationsUnderTest = new PipelineConfigurations();

    PipelineConfigurations.PipelineConfiguration missingDefaultDefinition =
        buildValidTestPipelineDefinition("array_imputation", 1);
    PipelineConfigurations.PipelineInputDefinitionConfig serviceProvidedInput =
        missingDefaultDefinition.getInputDefinitionConfigs().get(0);
    serviceProvidedInput.setUserProvided(false);
    serviceProvidedInput.setDefaultValue(null);

    pipelineConfigurationsUnderTest.setPipelines(
        Map.of(
            "arrayImputation",
            Map.of("1", missingDefaultDefinition),
            "lowPassImputation",
            Map.of("1", buildValidTestPipelineDefinition("low_pass_imputation", 1))));
    pipelineConfigurationsUnderTest.setPipelineQuotas(
        Map.of(
            "arrayImputation",
            buildValidTestQuotaConfig(),
            "lowPassImputation",
            buildValidTestQuotaConfig()));

    assertThrows(
        IllegalArgumentException.class,
        pipelineConfigurationsUnderTest::validateRuntimeConfiguration);
  }

  // TODO add test that runs validate, not just cast
  @Test
  void runtimeValidationWithWronglyTypedServiceProvidedDefaultThrows() {
    PipelineConfigurations pipelineConfigurationsUnderTest = new PipelineConfigurations();

    PipelineConfigurations.PipelineConfiguration wronglyTypedServiceProvidedInputConfig =
        buildValidTestPipelineDefinition("array_imputation", 1);
    PipelineConfigurations.PipelineInputDefinitionConfig serviceProvidedInput =
        wronglyTypedServiceProvidedInputConfig.getInputDefinitionConfigs().get(0);
    serviceProvidedInput.setUserProvided(false);
    serviceProvidedInput.setType(PipelineVariableTypesEnum.INTEGER);

    pipelineConfigurationsUnderTest.setPipelines(
        Map.of(
            "arrayImputation",
            Map.of("1", wronglyTypedServiceProvidedInputConfig),
            "lowPassImputation",
            Map.of("1", buildValidTestPipelineDefinition("low_pass_imputation", 1))));
    pipelineConfigurationsUnderTest.setPipelineQuotas(
        Map.of(
            "arrayImputation",
            buildValidTestQuotaConfig(),
            "lowPassImputation",
            buildValidTestQuotaConfig()));

    assertThrows(
        IllegalArgumentException.class,
        pipelineConfigurationsUnderTest::validateRuntimeConfiguration);
  }

  @Test
  void getPipelineConfigurationReturnsDefinition() {
    PipelineConfigurations.PipelineConfiguration definition =
        pipelineConfigurations.getPipelineConfiguration("array_imputation_v1");

    assertNotNull(definition);
    assertEquals("array_imputation", definition.getMetadata().getPipelineName());
    assertEquals(1, definition.getMetadata().getPipelineVersion());
  }

  @Test
  void getPipelineConfigurationThrowsForMissingDefinition() {
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelineConfigurations.getPipelineConfiguration("array_imputation_v999"));
  }

  @Test
  void allServiceProvidedInputsHaveDefaultValues() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputDefinitionConfigs().stream())
        .filter(input -> !input.getUserProvided())
        .forEach(input -> assertNotNull(input.getDefaultValue()));
  }

  @Test
  void allUserProvidedFileInputsHaveDefinedFileSuffixes() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputDefinitionConfigs().stream())
        .filter(input -> input.getUserProvided() && input.getType().isFileLike())
        .forEach(input -> assertNotNull(input.getFileSuffix()));
  }

  @Test
  void allDefaultValuesForPipelineInputsAreCorrectType() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputDefinitionConfigs().stream())
        .filter(input -> input.getDefaultValue() != null)
        .forEach(
            input -> {
              PipelineInputDefinition modelInputDefinition = toModelInputDefinition(input);
              assertNull(input.getType().validate(modelInputDefinition, input.getDefaultValue()));
              assertNotNull(
                  input
                      .getType()
                      .cast(input.getName(), input.getDefaultValue(), new TypeReference<>() {}));
            });
  }

  private Stream<PipelineConfigurations.PipelineConfiguration> allConfiguredPipelines() {
    return pipelineConfigurations.getPipelines().values().stream()
        .flatMap(v -> v.values().stream());
  }

  private PipelineInputDefinition toModelInputDefinition(
      PipelineConfigurations.PipelineInputDefinitionConfig inputDefinitionConfig) {
    return PipelineInputDefinition.builder()
        .name(inputDefinitionConfig.getName())
        .wdlVariableName(inputDefinitionConfig.getWdlVariableName())
        .displayName(inputDefinitionConfig.getDisplayName())
        .description(inputDefinitionConfig.getDescription())
        .type(inputDefinitionConfig.getType())
        .isRequired(inputDefinitionConfig.getIsRequired())
        .userProvided(inputDefinitionConfig.getUserProvided())
        .defaultValue(inputDefinitionConfig.getDefaultValue())
        .minValue(inputDefinitionConfig.getMinValue())
        .maxValue(inputDefinitionConfig.getMaxValue())
        .fileSuffix(inputDefinitionConfig.getFileSuffix())
        .build();
  }

  private PipelineConfigurations.PipelineConfiguration buildValidTestPipelineDefinition(
      String pipelineName, Integer pipelineVersion) {
    PipelineConfigurations.PipelineConfiguration definition =
        new PipelineConfigurations.PipelineConfiguration();

    PipelineConfigurations.PipelineMetadataConfig metadata =
        new PipelineConfigurations.PipelineMetadataConfig();
    metadata.setPipelineName(pipelineName);
    metadata.setPipelineVersion(pipelineVersion);
    metadata.setDisplayName("Display Name");
    metadata.setMemoryRetryMultiplier(BigDecimal.valueOf(1.0));
    definition.setMetadata(metadata);

    PipelineConfigurations.PipelineInputDefinitionConfig input =
        new PipelineConfigurations.PipelineInputDefinitionConfig();
    input.setName("input");
    input.setWdlVariableName("input");
    input.setType(bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING);
    input.setIsRequired(true);
    input.setUserProvided(true);
    input.setDefaultValue("validDefault");
    definition.setInputDefinitionConfigs(List.of(input));

    PipelineConfigurations.PipelineOutputDefinitionConfig output =
        new PipelineConfigurations.PipelineOutputDefinitionConfig();
    output.setName("output");
    output.setWdlVariableName("output");
    output.setType(bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING);
    output.setIsRequired(true);
    definition.setOutputDefinitionConfigs(List.of(output));

    return definition;
  }

  private PipelineConfigurations.PipelineQuotaConfig buildValidTestQuotaConfig() {
    PipelineConfigurations.PipelineQuotaConfig quota =
        new PipelineConfigurations.PipelineQuotaConfig();
    quota.setDefaultQuota(10);
    quota.setMinQuotaConsumed(1);
    quota.setQuotaUnits(bio.terra.pipelines.common.utils.QuotaUnitsEnum.SAMPLES);
    return quota;
  }

  // TESTING/VALIDATING REAL CONFIG VALUES
  @Nested
  @DisplayName("Service Config Validations")
  class ServicePipelinesConfigValidations {
    @Test
    void mainPipelinesConfigPassesRuntimeValidation() {
      PipelineConfigurations pipelineConfigurations = loadMainPipelineConfigurations();

      assertDoesNotThrow(pipelineConfigurations::validateRuntimeConfiguration);
    }

    @Test
    void mainPipelinesConfigPassesDetailedValidation() {
      PipelineConfigurations pipelineConfigurations = loadMainPipelineConfigurations();

      assertDoesNotThrow(pipelineConfigurations::validateDetailedConfiguration);
    }

    @Test
    void allDefinedDefaultValuesForPipelineInputsAreNonBlankInMainConfig() {
      PipelineConfigurations pipelineConfigurations = loadMainPipelineConfigurations();

      allConfiguredPipelines(pipelineConfigurations)
          .flatMap(config -> config.getInputDefinitionConfigs().stream())
          .filter(input -> input.getDefaultValue() != null)
          .forEach(
              input -> {
                String defaultValue = input.getDefaultValue();
                assertFalse(
                    defaultValue.isBlank(), "Default value is blank for " + input.getName());
              });
    }

    @Test
    void mainConfigContainsDefinitionsForAllKnownPipelines() {
      PipelineConfigurations pipelineConfigurations = loadMainPipelineConfigurations();
      Map<String, Map<String, PipelineConfigurations.PipelineConfiguration>> allPipelines =
          pipelineConfigurations.getPipelines();

      for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
        Map<String, PipelineConfigurations.PipelineConfiguration> versions =
            allPipelines.get(pipelineEnum.getConfigKeyValue());
        assertNotNull(versions, "Missing pipeline map for " + pipelineEnum.getConfigKeyValue());
        assertFalse(
            versions.isEmpty(), "No versions configured for " + pipelineEnum.getConfigKeyValue());
      }
    }

    private Stream<PipelineConfigurations.PipelineConfiguration> allConfiguredPipelines(
        PipelineConfigurations pipelineConfigurations) {
      return pipelineConfigurations.getPipelines().values().stream()
          .flatMap(v -> v.values().stream());
    }

    private PipelineConfigurations loadMainPipelineConfigurations() {
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
            .bind("pipelines.configurations", Bindable.of(PipelineConfigurations.class))
            .orElseThrow(
                () ->
                    new IllegalStateException(
                        "Could not bind pipelines.configurations from " + pipelinesConfigPath));
      } catch (IOException e) {
        throw new IllegalStateException("Failed to load main pipelines config", e);
      }
    }
  }
}
