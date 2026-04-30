package bio.terra.pipelines.configuration.internal;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.app.configuration.internal.PipelineCatalogConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelineCatalogConfigurationsTest extends BaseEmbeddedDbTest {

  @Autowired private PipelineCatalogConfigurations pipelineCatalogConfigurations;

  @Test
  void catalogConfigurationBindsFromYaml() {
    PipelineCatalogConfigurations.PipelineDefinition definition =
        pipelineCatalogConfigurations.getDefinitions().get("array_imputation").get("1");

    assertEquals("Test Array Imputation", definition.getDisplayName());
    assertEquals("imputation", definition.getPipelineType());
    assertEquals(2500, definition.getQuota().getDefaultQuota());
    assertEquals(QuotaUnitsEnum.SAMPLES, definition.getQuota().getQuotaUnits());
    assertEquals(1, definition.getInputDefinitions().size());
    assertEquals(1, definition.getOutputDefinitions().size());
  }

  @Test
  void validateWithDuplicateInputNamesThrows() {
    PipelineCatalogConfigurations config = new PipelineCatalogConfigurations();
    config.setDefinitions(
        Map.of("array_imputation", Map.of("1", buildDefinitionWithDuplicateInputNames())));

    assertThrows(IllegalArgumentException.class, config::validate);
  }

  @Test
  void validateWithInvalidInputRangeThrows() {
    PipelineCatalogConfigurations config = new PipelineCatalogConfigurations();
    config.setDefinitions(
        Map.of("array_imputation", Map.of("1", buildDefinitionWithInvalidInputRange())));

    assertThrows(IllegalArgumentException.class, config::validate);
  }

  @Test
  void validateWithValidDefinitionPasses() {
    PipelineCatalogConfigurations config = new PipelineCatalogConfigurations();
    config.setDefinitions(Map.of("array_imputation", Map.of("1", buildValidDefinition())));

    assertDoesNotThrow(config::validate);
  }

  private static PipelineCatalogConfigurations.PipelineDefinition buildValidDefinition() {
    PipelineCatalogConfigurations.PipelineDefinition definition =
        new PipelineCatalogConfigurations.PipelineDefinition();
    definition.setDisplayName("Array Imputation");
    definition.setDescription("Test description");
    definition.setPipelineType("imputation");
    definition.setToolName("ImputationBeagle");

    PipelineCatalogConfigurations.QuotaDefinition quota =
        new PipelineCatalogConfigurations.QuotaDefinition();
    quota.setDefaultQuota(2500);
    quota.setMinQuotaConsumed(500);
    quota.setQuotaUnits(QuotaUnitsEnum.SAMPLES);
    definition.setQuota(quota);

    PipelineCatalogConfigurations.InputDefinition input =
        new PipelineCatalogConfigurations.InputDefinition();
    input.setName("multiSampleVcf");
    input.setWdlVariableName("multi_sample_vcf");
    input.setDisplayName("multi-sample VCF");
    input.setDescription("test input");
    input.setType(PipelineVariableTypesEnum.FILE);
    input.setRequired(true);
    input.setUserProvided(true);
    input.setExpectsCustomValue(false);
    input.setFileSuffix(".vcf.gz");

    PipelineCatalogConfigurations.OutputDefinition output =
        new PipelineCatalogConfigurations.OutputDefinition();
    output.setName("imputedMultiSampleVcf");
    output.setWdlVariableName("imputed_multi_sample_vcf");
    output.setDisplayName("imputed VCF");
    output.setDescription("test output");
    output.setType(PipelineVariableTypesEnum.FILE);
    output.setRequired(true);

    definition.setInputDefinitions(List.of(input));
    definition.setOutputDefinitions(List.of(output));
    return definition;
  }

  private static PipelineCatalogConfigurations.PipelineDefinition
      buildDefinitionWithDuplicateInputNames() {
    PipelineCatalogConfigurations.PipelineDefinition definition = buildValidDefinition();

    PipelineCatalogConfigurations.InputDefinition duplicateInput =
        new PipelineCatalogConfigurations.InputDefinition();
    duplicateInput.setName("multiSampleVcf");
    duplicateInput.setWdlVariableName("some_other_wdl_name");
    duplicateInput.setDisplayName("duplicate input");
    duplicateInput.setDescription("duplicate input");
    duplicateInput.setType(PipelineVariableTypesEnum.FILE);
    duplicateInput.setRequired(true);
    duplicateInput.setUserProvided(true);
    duplicateInput.setExpectsCustomValue(false);
    duplicateInput.setFileSuffix(".vcf.gz");

    definition.setInputDefinitions(
        List.of(definition.getInputDefinitions().get(0), duplicateInput));
    return definition;
  }

  private static PipelineCatalogConfigurations.PipelineDefinition
      buildDefinitionWithInvalidInputRange() {
    PipelineCatalogConfigurations.PipelineDefinition definition = buildValidDefinition();

    PipelineCatalogConfigurations.InputDefinition rangeInput =
        new PipelineCatalogConfigurations.InputDefinition();
    rangeInput.setName("minDr2ForInclusion");
    rangeInput.setWdlVariableName("min_dr2_for_inclusion");
    rangeInput.setDisplayName("minimum dr2");
    rangeInput.setDescription("range input");
    rangeInput.setType(PipelineVariableTypesEnum.FLOAT);
    rangeInput.setRequired(false);
    rangeInput.setUserProvided(true);
    rangeInput.setExpectsCustomValue(false);
    rangeInput.setMinValue(1.0);
    rangeInput.setMaxValue(0.0);

    definition.setInputDefinitions(List.of(rangeInput));
    return definition;
  }
}
