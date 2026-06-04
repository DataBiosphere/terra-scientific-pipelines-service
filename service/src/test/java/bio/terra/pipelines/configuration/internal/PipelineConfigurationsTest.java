package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
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
  void testPipelinesCommonConfiguration() {
    PipelineConfigurations.CommonConfiguration commonConfiguration =
        pipelineConfigurations.getCommon();
    assertEquals(2, commonConfiguration.getUserDataTtlDays());

    assertEquals(1, commonConfiguration.getQuotaConsumedPollingIntervalSeconds());
    assertTrue(commonConfiguration.isQuotaConsumedUseCallCaching());

    assertEquals(1, commonConfiguration.getInputQcPollingIntervalSeconds());
    assertTrue(commonConfiguration.isInputQcUseCallCaching());

    assertEquals(2, commonConfiguration.getMainToolPollingIntervalSeconds());
    assertTrue(commonConfiguration.isMainToolUseCallCaching());
    assertFalse(commonConfiguration.isMainToolDeleteIntermediateFiles());

    assertEquals(
        "gs://test_bucket/test_path/to/monitoring/script.sh",
        commonConfiguration.getMonitoringScriptPath());
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
        100,
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

  @Test
  void allServiceProvidedInputsHaveDefaultValues() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputDefinitions().stream())
        .filter(input -> !input.getUserProvided())
        .forEach(input -> assertNotNull(input.getDefaultValue()));
  }

  @Test
  void allUserProvidedFileInputsHaveDefinedFileSuffixes() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputDefinitions().stream())
        .filter(input -> input.getUserProvided() && input.getType().isFileLike())
        .forEach(input -> assertNotNull(input.getFileSuffix()));
  }

  @Test
  void allDefaultValuesForPipelineInputsAreCorrectType() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputDefinitions().stream())
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

  private Stream<PipelineConfigurations.WdlBasedPipelineConfiguration> allConfiguredPipelines() {
    return pipelineConfigurations.getPipelines().values().stream()
        .flatMap(v -> v.values().stream());
  }

  private PipelineInputDefinition toModelInputDefinition(
      PipelineConfigurations.PipelineInputDefinitionConfiguration inputDefinitionConfig) {
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

  private PipelineConfigurations.WdlBasedPipelineConfiguration buildValidTestPipelineDefinition(
      String pipelineName, Integer pipelineVersion) {
    PipelineConfigurations.WdlBasedPipelineConfiguration definition =
        new PipelineConfigurations.WdlBasedPipelineConfiguration();
    definition.setDisplayName("Display Name");
    definition.setPipelineType("WDL");
    definition.setToolName("ToolName");
    definition.setMemoryRetryMultiplier(BigDecimal.valueOf(1.0));

    PipelineConfigurations.PipelineInputDefinitionConfiguration input =
        new PipelineConfigurations.PipelineInputDefinitionConfiguration();
    input.setName("input");
    input.setWdlVariableName("input");
    input.setType(bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING);
    input.setIsRequired(true);
    input.setUserProvided(true);
    input.setDefaultValue("validDefault");
    definition.setInputDefinitions(List.of(input));

    PipelineConfigurations.PipelineOutputDefinitionConfiguration output =
        new PipelineConfigurations.PipelineOutputDefinitionConfiguration();
    output.setName("output");
    output.setWdlVariableName("output");
    output.setType(bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING);
    output.setIsRequired(true);
    definition.setOutputDefinitions(List.of(output));

    return definition;
  }

  private PipelineConfigurations.PipelineQuotaConfiguration buildValidTestQuotaConfig() {
    PipelineConfigurations.PipelineQuotaConfiguration quota =
        new PipelineConfigurations.PipelineQuotaConfiguration();
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
    void servicePipelinesConfigPassesValidation() {
      PipelineConfigurations servicePipelineConfigurations = loadServicePipelineConfigurations();

      assertDoesNotThrow(servicePipelineConfigurations::validateConfiguration);
    }

    @Test
    void allDefinedDefaultValuesForPipelineInputsAreNonBlankInServiceConfig() {
      PipelineConfigurations servicePipelineConfigurations = loadServicePipelineConfigurations();

      allConfiguredPipelinesStream(servicePipelineConfigurations)
          .flatMap(config -> config.getInputDefinitions().stream())
          .filter(input -> input.getDefaultValue() != null)
          .forEach(
              input -> {
                String defaultValue = input.getDefaultValue();
                assertFalse(
                    defaultValue.isBlank(), "Default value is blank for " + input.getName());
              });
    }

    @Test
    void serviceConfigContainsDefinitionsForAllKnownPipelines() {
      PipelineConfigurations servicePipelineConfigurations = loadServicePipelineConfigurations();
      Map<String, Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration>> allPipelines =
          servicePipelineConfigurations.getPipelines();

      for (PipelinesEnum pipelineEnum : PipelinesEnum.values()) {
        Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versions =
            allPipelines.get(pipelineEnum.getConfigKeyValue());
        assertNotNull(versions, "Missing pipeline map for " + pipelineEnum.getConfigKeyValue());
        assertFalse(
            versions.isEmpty(), "No versions configured for " + pipelineEnum.getConfigKeyValue());
      }
    }

    private Stream<PipelineConfigurations.WdlBasedPipelineConfiguration>
        allConfiguredPipelinesStream(PipelineConfigurations pipelineConfigurations) {
      return pipelineConfigurations.getPipelines().values().stream()
          .flatMap(v -> v.values().stream());
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
  }
}
