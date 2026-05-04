package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.app.configuration.internal.PipelineCatalogConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.repositories.PipelineQuotasRepository;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelinesServiceDatabaseTest extends BaseEmbeddedDbTest {
  @Autowired PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;
  @Autowired PipelineQuotasRepository pipelineQuotasRepository;
  @Autowired PipelinesService pipelinesService;
  @Autowired PipelineCatalogService pipelineCatalogService;

  @Test
  void allPipelineEnumsExist() {
    // make sure all the pipelines in the enum exist in the table
    for (PipelinesEnum p : PipelinesEnum.values()) {
      assertTrue(pipelineRuntimeMetadataRepository.existsByNameAndHiddenIsFalse(p));
    }
  }

  @Test
  void allPipelinesHaveDefinedInputs() {
    // make sure all configured pipelines have defined inputs
    for (PipelinesEnum p : PipelinesEnum.values()) {
      Pipeline pipeline = pipelinesService.getPipeline(p, null, false);
      assertNotNull(pipeline.getPipelineInputDefinitions());
    }
  }

  @Test
  void allPipelineInputDefinitionsAreProperlyTyped() {
    // make sure all pipeline inputs have defined types matching the enum
    for (PipelineInputDefinition p : getAllPipelineInputDefinitions()) {
      assertDoesNotThrow(() -> p.getType());
    }
  }

  @Test
  void allDefaultValuesForPipelineInputsAreCorrectType() {
    // make sure all pipeline input definition default values pass type validation and are cast-able
    for (PipelineInputDefinition p : getAllPipelineInputDefinitions()) {
      if (p.getDefaultValue() != null) {
        PipelineVariableTypesEnum inputType = p.getType();
        assertNull(inputType.validate(p, p.getDefaultValue()));
        assertNotNull(inputType.cast(p.getName(), p.getDefaultValue(), new TypeReference<>() {}));
      }
    }
  }

  @Test
  void allOptionalInputsHaveDefaultValues() {
    // make sure all optional inputs have default values
    for (PipelineInputDefinition p : getAllPipelineInputDefinitions()) {
      if (!p.isRequired()) {
        assertNotNull(p.getDefaultValue());
      }
    }
  }

  @Test
  void allServiceProvidedInputsWithoutCustomValuesHaveDefaultValues() {
    // make sure all service-provided inputs that aren't marked as having custom values, have
    // default values
    for (PipelineInputDefinition p : getAllPipelineInputDefinitions()) {
      if (!p.isUserProvided() && !p.isExpectsCustomValue()) {
        assertNotNull(p.getDefaultValue());
      }
    }
  }

  @Test
  void allFileInputsHaveDefinedFileSuffixes() {
    // make sure each FILE input has a defined value for file_suffix
    for (PipelineInputDefinition p : getAllPipelineInputDefinitions()) {
      if (p.getType() == PipelineVariableTypesEnum.FILE
          || p.getType() == PipelineVariableTypesEnum.FILE_ARRAY) {
        assertNotNull(p.getFileSuffix());
      }
    }
  }

  @Test
  void imputationPipelineV1HasCorrectInputsAndOutputs() {
    Pipeline pipeline = pipelinesService.getPipeline(PipelinesEnum.ARRAY_IMPUTATION, 1, true);
    assertPipelineMatchesCatalogDefinitions(pipeline, PipelinesEnum.ARRAY_IMPUTATION, 1);
  }

  @Test
  void imputationPipelineV2HasCorrectInputsAndOutputs() {
    Pipeline pipeline = pipelinesService.getPipeline(PipelinesEnum.ARRAY_IMPUTATION, 2, true);
    assertPipelineMatchesCatalogDefinitions(pipeline, PipelinesEnum.ARRAY_IMPUTATION, 2);
  }

  @Test
  void allPipelinesHaveDefinedOutputs() {
    // make sure all configured pipelines have defined outputs
    for (PipelinesEnum p : PipelinesEnum.values()) {
      Pipeline pipeline = pipelinesService.getPipeline(p, null, false);
      assertNotNull(pipeline.getPipelineOutputDefinitions());
    }
  }

  @Test
  void imputationPipelineHasCorrectOutputs() {
    Pipeline pipeline = pipelinesService.getPipeline(PipelinesEnum.ARRAY_IMPUTATION, 1, false);
    assertPipelineMatchesCatalogDefinitions(pipeline, PipelinesEnum.ARRAY_IMPUTATION, 1);
  }

  @Test
  void allPipelinesInPipelineQuotasHaveDefinedPipelines() {
    // make sure all the pipelines in the pipeline quotas table have defined pipelines in the
    // pipelines table
    for (PipelineQuota pq : pipelineQuotasRepository.findAll()) {
      assertTrue(
          pipelineRuntimeMetadataRepository.existsByNameAndHiddenIsFalse(pq.getPipelineName()));
    }
  }

  private List<PipelineInputDefinition> getAllPipelineInputDefinitions() {
    return pipelinesService.getPipelines(true).stream()
        .flatMap(pipeline -> pipeline.getPipelineInputDefinitions().stream())
        .toList();
  }

  private void assertPipelineMatchesCatalogDefinitions(
      Pipeline pipeline, PipelinesEnum pipelineName, Integer version) {
    PipelineCatalogConfigurations.PipelineDefinition definition =
        pipelineCatalogService
            .getDefinition(pipelineName, version)
            .orElseThrow(() -> new AssertionError("Missing catalog definition for pipeline test"));

    assertEquals(
        definition.getInputDefinitions().size(), pipeline.getPipelineInputDefinitions().size());
    assertEquals(
        definition.getInputDefinitions().stream()
            .map(PipelineCatalogConfigurations.VariableDefinition::getName)
            .collect(Collectors.toSet()),
        pipeline.getPipelineInputDefinitions().stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet()));
    assertEquals(
        definition.getInputDefinitions().stream()
            .map(PipelineCatalogConfigurations.VariableDefinition::getWdlVariableName)
            .collect(Collectors.toSet()),
        pipeline.getPipelineInputDefinitions().stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet()));

    assertEquals(
        definition.getOutputDefinitions().size(), pipeline.getPipelineOutputDefinitions().size());
    assertEquals(
        definition.getOutputDefinitions().stream()
            .map(PipelineCatalogConfigurations.VariableDefinition::getName)
            .collect(Collectors.toSet()),
        pipeline.getPipelineOutputDefinitions().stream()
            .map(PipelineOutputDefinition::getName)
            .collect(Collectors.toSet()));
    assertEquals(
        definition.getOutputDefinitions().stream()
            .map(PipelineCatalogConfigurations.VariableDefinition::getWdlVariableName)
            .collect(Collectors.toSet()),
        pipeline.getPipelineOutputDefinitions().stream()
            .map(PipelineOutputDefinition::getWdlVariableName)
            .collect(Collectors.toSet()));

    assertEquals(
        Set.of(pipeline.getId()),
        pipeline.getPipelineInputDefinitions().stream()
            .map(PipelineInputDefinition::getPipelineId)
            .collect(Collectors.toSet()));
    assertEquals(
        Set.of(pipeline.getId()),
        pipeline.getPipelineOutputDefinitions().stream()
            .map(PipelineOutputDefinition::getPipelineId)
            .collect(Collectors.toSet()));
  }
}
