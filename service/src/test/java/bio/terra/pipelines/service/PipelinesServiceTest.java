package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelineInputDefinitionsRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.util.*;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelinesServiceTest extends BaseEmbeddedDbTest {
  @Autowired PipelinesService pipelinesService;
  @Autowired PipelinesRepository pipelinesRepository;
  @Autowired PipelineInputDefinitionsRepository pipelineInputDefinitionsRepository;

  @Test
  void getCorrectNumberOfPipelines() {
    // migrations insert one pipeline (imputation) so make sure we find it
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    UUID workspaceId = UUID.randomUUID();

    pipelinesRepository.save(
        new Pipeline(
            "pipelineName",
            "1.0.0",
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "wdlMethodName",
            workspaceId,
            null));

    pipelineList = pipelinesService.getPipelines();
    assertEquals(2, pipelineList.size());
    Pipeline savedPipeline = pipelineList.get(1);
    assertEquals("pipelineName", savedPipeline.getName());
    assertEquals("1.0.0", savedPipeline.getVersion());
    assertEquals("pipelineDisplayName", savedPipeline.getDisplayName());
    assertEquals("description", savedPipeline.getDescription());
    assertEquals("pipelineType", savedPipeline.getPipelineType());
    assertEquals("wdlUrl", savedPipeline.getWdlUrl());
    assertEquals("wdlMethodName", savedPipeline.getWdlMethodName());
    assertEquals(workspaceId, savedPipeline.getWorkspaceId());
  }

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
  void imputationPipelineHasCorrectInputs() {
    // make sure the imputation pipeline has the correct inputs
    Pipeline pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());

    List<PipelineInputDefinition> pipelineInputDefinitions = pipeline.getPipelineInputDefinitions();

    // currently we have one input for the imputation pipeline
    assertEquals(1, pipelineInputDefinitions.size());

    PipelineInputDefinition input1 = pipelineInputDefinitions.get(0);
    assertEquals("multi_sample_vcf", input1.getName());
    assertEquals("String", input1.getType());
    assertTrue(input1.getIsRequired());
    assertNotNull(input1.getId());
    // make sure the inputs are associated with the correct pipeline
    assertEquals(pipeline.getId(), input1.getPipelineId());
  }

  @Test
  void addPipelineInput() {
    Pipeline pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());
    List<PipelineInputDefinition> pipelineInputDefinitions = pipeline.getPipelineInputDefinitions();
    assertEquals(1, pipelineInputDefinitions.size());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput = new PipelineInputDefinition();
    newInput.setPipelineId(pipeline.getId());
    newInput.setName("newInput");
    newInput.setType("Int");
    newInput.setIsRequired(false);

    pipelineInputDefinitionsRepository.save(newInput);

    pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());
    pipelineInputDefinitions = pipeline.getPipelineInputDefinitions();
    assertEquals(2, pipelineInputDefinitions.size());

    PipelineInputDefinition savedInput = pipelineInputDefinitions.get(1);
    assertEquals("newInput", savedInput.getName());
    assertEquals("Int", savedInput.getType());
    assertFalse(savedInput.getIsRequired());
    assertEquals(pipeline.getId(), savedInput.getPipelineId());
  }

  @Test
  void testPipelineToString() {
    // test .ToString() method on Pipeline Entity
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      assertEquals(
          String.format(
              "Pipeline[pipelineName=%s, version=%s, displayName=%s, description=%s, pipelineType=%s, wdlUrl=%s, wdlMethodName=%s, workspaceId=%s]",
              p.getName(),
              p.getVersion(),
              p.getDisplayName(),
              p.getDescription(),
              p.getPipelineType(),
              p.getWdlUrl(),
              p.getWdlMethodName(),
              p.getWorkspaceId()),
          p.toString());
    }
  }

  @Test
  void testPipelineHashCode() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      assertEquals(
          new HashCodeBuilder(17, 31)
              .append(p.getId())
              .append(p.getName())
              .append(p.getVersion())
              .append(p.getDisplayName())
              .append(p.getDescription())
              .append(p.getPipelineType())
              .append(p.getWdlUrl())
              .append(p.getWdlMethodName())
              .append(p.getWorkspaceId())
              .toHashCode(),
          p.hashCode());
    }
  }

  @Test
  void updatePipelineWorkspaceId() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum);
    UUID savedWorkspaceId = UUID.randomUUID();
    // make sure the current pipeline does not have the workspace id we're trying to update with
    assertNotEquals(savedWorkspaceId, p.getWorkspaceId());

    // update pipeline workspace id
    pipelinesService.updatePipelineWorkspaceId(pipelinesEnum, savedWorkspaceId);
    p = pipelinesService.getPipeline(pipelinesEnum);

    // assert the workspace id has been updated
    assertEquals(savedWorkspaceId, p.getWorkspaceId());
  }

  // input validation tests

  @Test
  void validateInputsWithStringAndInt() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of("multi_sample_vcf", "this is a string", "new_integer_input", 123));
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsWithStringAndStringyInt() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of("multi_sample_vcf", "this is a string", "new_integer_input", "123"));
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsEmpty() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Object emptyInputs = new Object();
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, emptyInputs));
  }

  @Test
  void validateInputsMultipleEmpty() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object emptyInputs = new Object();
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, emptyInputs));
  }

  @Test
  void validateInputsNotAMap() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Object inputs = new ArrayList<>(List.of("this", "is", "not", "a", "string"));
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsNotAString() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of("multi_sample_vcf", List.of("this", "is", "not", "a", "string")));
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsArrayString() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(
            pipeline.getId(), "new_array_string_input", "Array_String", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of(
                "multi_sample_vcf",
                "this is a string",
                "new_array_string_input",
                List.of("this", "is", "an", "array", "of", "strings")));
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsArrayNotAnArray() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(
            pipeline.getId(), "new_array_string_input", "Array_String", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of(
                "multi_sample_vcf",
                "this is a string",
                "new_array_string_input",
                "this is not an array"));
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsStringNotAnInteger() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of(
                "multi_sample_vcf",
                "this is a string",
                "new_integer_input",
                "this is not an integer"));
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsFloatNotAnInteger() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of("multi_sample_vcf", "this is a string", "new_integer_input", "2.3"));
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateOptionalInputIsValidType() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", false);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of("multi_sample_vcf", "this is a string", "new_integer_input", 123));
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateOptionalInputBadType() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", false);

    pipelineInputDefinitionsRepository.save(newInput);

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of(
                "multi_sample_vcf",
                "this is a string",
                "new_integer_input",
                "this is not an integer"));
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsNoneRequired() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;

    // remove all pipeline inputs definitions
    pipelineInputDefinitionsRepository.deleteAll();

    Object inputs =
        new LinkedHashMap<String, Object>(
            Map.of(
                "multi_sample_vcf",
                "this is a string",
                "new_integer_input",
                "this is not an integer"));
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));

    Object inputsEmpty = new Object();
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputsEmpty));
  }
}
