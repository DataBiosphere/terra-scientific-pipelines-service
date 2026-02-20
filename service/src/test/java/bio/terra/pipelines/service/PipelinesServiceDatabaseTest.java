package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.repositories.PipelineInputDefinitionsRepository;
import bio.terra.pipelines.db.repositories.PipelineQuotasRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;

class PipelinesServiceDatabaseTest extends BaseEmbeddedDbTest {
  @Autowired PipelinesRepository pipelinesRepository;
  @Autowired PipelineInputDefinitionsRepository pipelineInputDefinitionsRepository;
  @Autowired PipelineQuotasRepository pipelineQuotasRepository;

  @Test
  void allPipelineEnumsExist() {
    // make sure all the pipelines in the enum exist in the table
    for (PipelinesEnum p : PipelinesEnum.values()) {
      assertTrue(pipelinesRepository.existsByNameAndHiddenIsFalse(p));
    }
  }

  @Test
  void allPipelinesHaveDefinedInputs() {
    // make sure all the pipelines in the enum have defined inputs
    for (PipelinesEnum p : PipelinesEnum.values()) {
      Pipeline pipeline = pipelinesRepository.findByNameAndHiddenIsFalse(p);
      assertNotNull(pipeline.getPipelineInputDefinitions());
    }
  }

  @Test
  void allPipelineInputDefinitionsAreProperlyTyped() {
    // make sure all pipeline inputs have defined types matching the enum
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      assertDoesNotThrow(() -> p.getType());
    }
  }

  @Test
  void allDefaultValuesForPipelineInputsAreCorrectType() {
    // make sure all pipeline input definition default values pass type validation and are cast-able
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      if (p.getDefaultValue() != null) {
        PipelineVariableTypesEnum inputType = p.getType();
        assertNull(inputType.validate(p, p.getDefaultValue()));
        assertNotNull(inputType.cast(p.getName(), p.getDefaultValue(), new TypeReference<>() {}));
      }
    }
  }

  @Test
  void allServiceProvidedInputsWithoutCustomValuesHaveDefaultValues() {
    // make sure all service-provided inputs that aren't marked as having custom values, have
    // default values
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      if (!p.isUserProvided() && !p.isExpectsCustomValue()) {
        assertNotNull(p.getDefaultValue());
      }
    }
  }

  @Test
  void allFileInputsHaveDefinedFileSuffixes() {
    // make sure each FILE input has a defined value for file_suffix
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      if (p.getType() == PipelineVariableTypesEnum.FILE
          || p.getType() == PipelineVariableTypesEnum.FILE_ARRAY) {
        assertNotNull(p.getFileSuffix());
      }
    }
  }

  @Test
  void imputationPipelineHasCorrectInputs() {
    // make sure the imputation pipeline has the correct inputs
    Pipeline pipeline =
        pipelinesRepository.findByNameAndHiddenIsFalse(PipelinesEnum.ARRAY_IMPUTATION);

    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipeline.getPipelineInputDefinitions();

    // there should be 3 user-provided inputs and 5 service-provided inputs
    assertEquals(
        3,
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .count());
    assertEquals(
        5,
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::isUserProvided))
            .count());

    // check user-provided inputs
    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .toList()
            .stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(Set.of("multi_sample_vcf", "output_basename", "min_dr2_for_inclusion")));

    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .toList()
            .stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(Set.of("multiSampleVcf", "outputBasename", "minDr2ForInclusion")));

    // check service-provided inputs
    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::isUserProvided))
            .toList()
            .stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of(
                    "contigs",
                    "genetic_maps_path",
                    "ref_dict",
                    "reference_panel_path_prefix",
                    "pipeline_header_line")));

    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::isUserProvided))
            .toList()
            .stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of(
                    "contigs",
                    "geneticMapsPath",
                    "refDict",
                    "referencePanelPathPrefix",
                    "pipelineHeaderLine")));

    // make sure the inputs are associated with the correct pipeline
    assertEquals(
        Set.of(pipeline.getId()),
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getPipelineId)
            .collect(Collectors.toSet()));
  }

  @Test
  void addDuplicatePipelineInputThrows() {
    Pipeline pipeline =
        pipelinesRepository.findByNameAndHiddenIsFalse(PipelinesEnum.ARRAY_IMPUTATION);

    // add a pipeline input that already exists
    PipelineInputDefinition newInput = new PipelineInputDefinition();
    newInput.setPipelineId(pipeline.getId());
    newInput.setName("multiSampleVcf");
    newInput.setWdlVariableName("multi_sample_vcf");
    newInput.setType(PipelineVariableTypesEnum.INTEGER);
    newInput.setRequired(false);
    newInput.setUserProvided(true);
    newInput.setDefaultValue("42");

    assertThrows(
        DataIntegrityViolationException.class,
        () -> pipelineInputDefinitionsRepository.save(newInput));
  }

  @Test
  void allPipelinesHaveDefinedOutputs() {
    // make sure all the pipelines in the enum have defined outputs
    for (PipelinesEnum p : PipelinesEnum.values()) {
      Pipeline pipeline = pipelinesRepository.findByNameAndHiddenIsFalse(p);
      assertNotNull(pipeline.getPipelineOutputDefinitions());
    }
  }

  @Test
  void imputationPipelineHasCorrectOutputs() {
    // make sure the imputation pipeline has the correct outputs
    Pipeline pipeline =
        pipelinesRepository.findByNameAndHiddenIsFalse(PipelinesEnum.ARRAY_IMPUTATION);

    List<PipelineOutputDefinition> allPipelineOutputDefinitions =
        pipeline.getPipelineOutputDefinitions();

    // there should be 4 outputs
    assertEquals(4, allPipelineOutputDefinitions.stream().count());

    // check outputs
    assertTrue(
        allPipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of(
                    "imputed_multi_sample_vcf",
                    "imputed_multi_sample_vcf_index",
                    "chunks_info",
                    "contigs_info")));

    assertTrue(
        allPipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of(
                    "imputedMultiSampleVcf",
                    "imputedMultiSampleVcfIndex",
                    "chunksInfo",
                    "contigsInfo")));

    // make sure the outputs are associated with the correct pipeline
    assertEquals(
        Set.of(pipeline.getId()),
        allPipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getPipelineId)
            .collect(Collectors.toSet()));
  }

  @Test
  void allPipelinesInPipelineQuotasHaveDefinedPipelines() {
    // make sure all the pipelines in the pipeline quotas table have defined pipelines in the
    // pipelines table
    for (PipelineQuota pq : pipelineQuotasRepository.findAll()) {
      assertTrue(pipelinesRepository.existsByNameAndHiddenIsFalse(pq.getPipelineName()));
    }
  }
}
