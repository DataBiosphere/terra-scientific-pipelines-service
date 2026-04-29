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
    // make sure all the pipelines in the enum exist in the database with at least one version
    for (PipelinesEnum p : PipelinesEnum.values()) {
      assertNotNull(pipelinesRepository.findAllByNameOrderByVersionDesc(p));
    }
  }

  @Test
  void allPipelinesHaveDefinedInputs() {
    // make sure all the pipelines in the enum have defined inputs
    for (PipelinesEnum p : PipelinesEnum.values()) {
      List<Pipeline> pipelines = pipelinesRepository.findAllByNameOrderByVersionDesc(p);
      for (Pipeline pipeline : pipelines) {
        assertNotNull(pipeline.getPipelineInputDefinitions());
      }
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
  void allServiceProvidedInputsWithCustomValuesAreRequired() {
    // make sure all service-provided inputs that are marked as having custom values, are required.
    // otherwise if they are marked as not required and have no default value, the logic in
    // PipelineInputsOutputsService.formatPipelineInputs() will skip the input and not do the proper
    // substitution.
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      if (!p.isUserProvided() && p.isExpectsCustomValue()) {
        assertTrue(p.isRequired());
      }
    }
  }

  @Test
  void allFileInputsHaveDefinedFileSuffixes() {
    // make sure each user-provided FILE input has a defined value for file_suffix
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      if (p.isUserProvided() && p.getType().isFileLike()) {
        assertNotNull(p.getFileSuffix());
      }
    }
  }

  @Test
  void arrayImputationPipelineV1HasCorrectInputsAndOutputs() {
    // make sure the imputation pipeline has the correct inputs
    Pipeline pipeline = pipelinesRepository.findByNameAndVersion(PipelinesEnum.ARRAY_IMPUTATION, 1);

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

    List<PipelineOutputDefinition> pipelineOutputDefinitions =
        pipeline.getPipelineOutputDefinitions();
    assertEquals(4, pipelineOutputDefinitions.size());
    assertEquals(
        Set.of("imputedMultiSampleVcf", "imputedMultiSampleVcfIndex", "contigsInfo", "chunksInfo"),
        pipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getName)
            .collect(Collectors.toSet()));
  }

  @Test
  void arrayImputationPipelineV2HasCorrectInputsAndOutputs() {
    // make sure the imputation pipeline has the correct inputs
    Pipeline pipeline = pipelinesRepository.findByNameAndVersion(PipelinesEnum.ARRAY_IMPUTATION, 2);

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

    List<PipelineOutputDefinition> pipelineOutputDefinitions =
        pipeline.getPipelineOutputDefinitions();
    assertEquals(6, pipelineOutputDefinitions.size());
    assertEquals(
        Set.of(
            "imputedMultiSampleVcf",
            "imputedMultiSampleVcfIndex",
            "imputedHomRefSitesOnlyVcf",
            "imputedHomRefSitesOnlyVcfIndex",
            "contigsInfo",
            "chunksInfo"),
        pipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getName)
            .collect(Collectors.toSet()));
  }

  @Test
  void lowPassImputationPipelineV1HasCorrectInputsAndOutputs() {
    Pipeline pipeline =
        pipelinesRepository.findByNameAndVersion(PipelinesEnum.LOW_PASS_IMPUTATION, 1);

    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipeline.getPipelineInputDefinitions();

    // there should be 5 user-provided inputs and 5 service-provided inputs
    assertEquals(
        5,
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
            .containsAll(
                Set.of("crams", "cram_indices", "sample_ids", "cram_manifest", "output_basename")));

    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .toList()
            .stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of("crams", "cramIndices", "sampleIds", "cramManifest", "outputBasename")));

    // check service-provided inputs
    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::isUserProvided))
            .toList()
            .stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of("contigs", "ref_dict", "reference_panel_prefix", "fasta", "fasta_index")));

    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::isUserProvided))
            .toList()
            .stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of("contigs", "refDict", "referencePanelPrefix", "fasta", "fastaIndex")));

    // make sure the inputs are associated with the correct pipeline
    assertEquals(
        Set.of(pipeline.getId()),
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getPipelineId)
            .collect(Collectors.toSet()));

    // check outputs
    List<PipelineOutputDefinition> allPipelineOutputDefinitions =
        pipeline.getPipelineOutputDefinitions();

    // there should be 5 outputs
    assertEquals(5, allPipelineOutputDefinitions.stream().count());
    assertTrue(
        allPipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of(
                    "imputed_vcf",
                    "imputed_vcf_index",
                    "imputed_hom_ref_sites_only_vcf",
                    "imputed_hom_ref_sites_only_vcf_index",
                    "qc_metrics")));

    assertTrue(
        allPipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of(
                    "imputedVcf",
                    "imputedVcfIndex",
                    "imputedHomRefSitesOnlyVcf",
                    "imputedHomRefSitesOnlyVcfIndex",
                    "qcMetrics")));

    // make sure the outputs are associated with the correct pipeline
    assertEquals(
        Set.of(pipeline.getId()),
        allPipelineOutputDefinitions.stream()
            .map(PipelineOutputDefinition::getPipelineId)
            .collect(Collectors.toSet()));
  }

  @Test
  void addDuplicatePipelineInputThrows() {
    Pipeline pipeline = pipelinesRepository.findByNameAndVersion(PipelinesEnum.ARRAY_IMPUTATION, 1);

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
    // make sure all the pipelines in the enum have defined outputs for each version
    for (PipelinesEnum p : PipelinesEnum.values()) {
      List<Pipeline> pipelines = pipelinesRepository.findAllByOrderByNameAscVersionDesc();
      for (Pipeline pipeline : pipelines) {
        if (pipeline.getName() == p) {
          assertNotNull(pipeline.getPipelineOutputDefinitions());
        }
      }
    }
  }

  @Test
  void allPipelinesInPipelineQuotasHaveDefinedPipelines() {
    // make sure all the pipelines in the pipeline quotas table have defined pipelines in the
    // pipelines table
    for (PipelineQuota pq : pipelineQuotasRepository.findAll()) {
      assertNotNull(pipelinesRepository.findAllByNameOrderByVersionDesc(pq.getPipelineName()));
    }
  }
}
