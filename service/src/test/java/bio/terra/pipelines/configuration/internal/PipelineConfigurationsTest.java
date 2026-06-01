package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.fasterxml.jackson.core.type.TypeReference;
import java.math.BigDecimal;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelineConfigurationsTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineConfigurations pipelineConfigurations;

  @Test
  void testArrayImputationV1Configuration() {
    // note these are the values in test/resources/pipelines-config.yml and not production values
    PipelineConfigurations.PipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getArrayImputation().get("1");
    PipelineConfigurations.PipelineMetadataConfig metadata = pipelineConfiguration.getMetadata();

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(0.0);

    assertEquals(memoryRetryMultiplier, metadata.getMemoryRetryMultiplier());
  }

  @Test
  void testArrayImputationV2Configuration() {
    // note these are the values in test/resources/pipelines-config.yml and not production values
    PipelineConfigurations.PipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getArrayImputation().get("2");
    PipelineConfigurations.PipelineMetadataConfig metadata = pipelineConfiguration.getMetadata();

    BigDecimal memoryRetryMultiplier = BigDecimal.valueOf(1.4);

    assertEquals(memoryRetryMultiplier, metadata.getMemoryRetryMultiplier());
  }

  @Test
  void testLowPassImputationV1Configuration() {
    // note these are the values in test/resources/pipelines-config.yml and not production values
    PipelineConfigurations.PipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getLowPassImputation().get("1");
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
        pipelineConfigurations.getArrayImputation().get("1");

    assertNotNull(pipelineDefinitionConfig);
    assertEquals("array_imputation", pipelineDefinitionConfig.getMetadata().getPipelineName());
    assertEquals(1, pipelineDefinitionConfig.getMetadata().getPipelineVersion());
    assertEquals(1, pipelineDefinitionConfig.getInputs().size());
    assertEquals(1, pipelineDefinitionConfig.getOutputs().size());
    assertEquals(
        100,
        pipelineConfigurations
            .getQuotaForPipeline(PipelinesEnum.ARRAY_IMPUTATION)
            .getDefaultQuota());
  }

  @Test
  void invalidPipelineKeyInPipelineDefinitionsThrows() {
    PipelineConfigurations pipelineConfigurationsUnderTest = new PipelineConfigurations();
    pipelineConfigurationsUnderTest.setArrayImputation(
        Map.of("BadVersion", buildValidTestPipelineDefinition("array_imputation", 1)));

    assertThrows(
        IllegalArgumentException.class, pipelineConfigurationsUnderTest::validateConfiguration);
  }

  @Test
  void missingPipelineQuotaThrows() {
    PipelineConfigurations pipelineConfigurationsUnderTest = new PipelineConfigurations();
    pipelineConfigurationsUnderTest.setArrayImputation(
        Map.of("1", buildValidTestPipelineDefinition("array_imputation", 1)));
    pipelineConfigurationsUnderTest.setLowPassImputation(
        Map.of("1", buildValidTestPipelineDefinition("low_pass_imputation", 1)));
    pipelineConfigurationsUnderTest.setPipelineQuotas(
        Map.of("array_imputation", buildValidTestQuotaConfig()));

    assertThrows(
        IllegalArgumentException.class, pipelineConfigurationsUnderTest::validateConfiguration);
  }

  @Test
  void getWdlBasedPipelineConfigByKeyReturnsDefinition() {
    PipelineConfigurations.PipelineConfiguration definition =
        pipelineConfigurations.getWdlBasedPipelineConfigByKey("array_imputation_v1");

    assertNotNull(definition);
    assertEquals("array_imputation", definition.getMetadata().getPipelineName());
    assertEquals(1, definition.getMetadata().getPipelineVersion());
  }

  @Test
  void getWdlBasedPipelineConfigByKeyThrowsForMissingDefinition() {
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelineConfigurations.getWdlBasedPipelineConfigByKey("array_imputation_v999"));
  }

  @Test
  void allServiceProvidedInputsWithoutCustomValuesHaveDefaultValues() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputs().stream())
        .filter(input -> !input.getUserProvided() && !input.getExpectsCustomValue())
        .forEach(input -> assertNotNull(input.getDefaultValue()));
  }

  @Test
  void allServiceProvidedInputsWithCustomValuesAreRequired() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputs().stream())
        .filter(input -> !input.getUserProvided() && input.getExpectsCustomValue())
        .forEach(input -> assertTrue(input.getIsRequired()));
  }

  @Test
  void allUserProvidedFileInputsHaveDefinedFileSuffixes() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputs().stream())
        .filter(input -> input.getUserProvided() && input.getType().isFileLike())
        .forEach(input -> assertNotNull(input.getFileSuffix()));
  }

  @Test
  void allDefaultValuesForPipelineInputsAreCorrectType() {
    allConfiguredPipelines()
        .flatMap(config -> config.getInputs().stream())
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
    return Stream.concat(
        pipelineConfigurations.getArrayImputation().values().stream(),
        pipelineConfigurations.getLowPassImputation().values().stream());
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
    input.setExpectsCustomValue(false);
    definition.setInputs(List.of(input));

    PipelineConfigurations.PipelineOutputDefinitionConfig output =
        new PipelineConfigurations.PipelineOutputDefinitionConfig();
    output.setName("output");
    output.setWdlVariableName("output");
    output.setType(bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING);
    output.setIsRequired(true);
    definition.setOutputs(List.of(output));

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
}
