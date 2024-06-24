package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelineInputDefinitionsRepository;
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

  @Test
  void allPipelineEnumsExist() {
    // make sure all the pipelines in the enum exist in the table
    for (PipelinesEnum p : PipelinesEnum.values()) {
      assertTrue(pipelinesRepository.existsByName(p.getValue()));
    }
  }

  @Test
  void allPipelinesHaveDefinedInputs() {
    // make sure all the pipelines in the enum have defined inputs
    for (PipelinesEnum p : PipelinesEnum.values()) {
      Pipeline pipeline = pipelinesRepository.findByName(p.getValue());
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
        PipelineInputTypesEnum inputType = p.getType();
        assertNull(inputType.validate(p.getName(), p.getDefaultValue()));
        assertNotNull(inputType.cast(p.getName(), p.getDefaultValue(), new TypeReference<>() {}));
      }
    }
  }

  @Test
  void allOptionalAndServiceProvidedInputsHaveDefaultValues() {
    // make sure all optional and service-provided inputs have default values
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      if (!p.getIsRequired() || !p.getUserProvided()) {
        assertNotNull(p.getDefaultValue());
      }
    }
  }

  @Test
  void imputationPipelineHasCorrectInputs() {
    // make sure the imputation pipeline has the correct inputs
    Pipeline pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());

    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipeline.getPipelineInputDefinitions();

    // there should be 2 user-provided inputs and 4 service-provided inputs
    assertEquals(
        2,
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .count());
    assertEquals(
        4,
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::getUserProvided))
            .count());

    // check user-provided inputs
    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .toList()
            .stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(Set.of("multi_sample_vcf", "output_basename")));

    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .toList()
            .stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(Set.of("multiSampleVcf", "outputBasename")));

    // check service-provided inputs
    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::getUserProvided))
            .toList()
            .stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of("contigs", "genetic_maps_path", "ref_dict", "reference_panel_path_prefix")));

    assertTrue(
        allPipelineInputDefinitions.stream()
            .filter(Predicate.not(PipelineInputDefinition::getUserProvided))
            .toList()
            .stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())
            .containsAll(
                Set.of("contigs", "geneticMapsPath", "refDict", "referencePanelPathPrefix")));

    // make sure the inputs are associated with the correct pipeline
    assertEquals(
        Set.of(pipeline.getId()),
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getPipelineId)
            .collect(Collectors.toSet()));
  }

  @Test
  void addDuplicatePipelineInputThrows() {
    Pipeline pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());

    // add a pipeline input that already exists
    PipelineInputDefinition newInput = new PipelineInputDefinition();
    newInput.setPipelineId(pipeline.getId());
    newInput.setName("multiSampleVcf");
    newInput.setWdlVariableName("multi_sample_vcf");
    newInput.setType(PipelineInputTypesEnum.INTEGER);
    newInput.setIsRequired(false);
    newInput.setUserProvided(true);
    newInput.setDefaultValue("42");

    assertThrows(
        DataIntegrityViolationException.class,
        () -> pipelineInputDefinitionsRepository.save(newInput));
  }
}
