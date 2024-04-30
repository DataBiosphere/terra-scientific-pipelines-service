package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelineInputDefinitionsRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.util.*;
import java.util.stream.Stream;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.runner.RunWith;
import org.mockito.junit.MockitoJUnitRunner;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;

@RunWith(MockitoJUnitRunner.class)
class PipelinesServiceTest extends BaseEmbeddedDbTest {
  private Logger loggerMock;
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
  void allPipelineInputsAreProperlyTyped() {
    // make sure all pipeline inputs have defined types matching the enum
    for (PipelineInputDefinition p : pipelineInputDefinitionsRepository.findAll()) {
      assertDoesNotThrow(() -> PipelineInputTypesEnum.valueOf(p.getType().toUpperCase()));
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
    assertEquals("VCF", input1.getType());
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
    newInput.setType("Integer");
    newInput.setIsRequired(false);

    pipelineInputDefinitionsRepository.save(newInput);

    pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());
    pipelineInputDefinitions = pipeline.getPipelineInputDefinitions();
    assertEquals(2, pipelineInputDefinitions.size());

    PipelineInputDefinition savedInput = pipelineInputDefinitions.get(1);
    assertEquals("newInput", savedInput.getName());
    assertEquals("Integer", savedInput.getType());
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
            Map.of("multi_sample_vcf", "this/is/a/vcf/path.vcf.gz", "new_integer_input", 123));
    assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateInputsEmpty() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;

    Object emptyInputs = new Object();
    ValidationException exception =
        assertThrows(
            ValidationException.class,
            () -> pipelinesService.validateInputs(pipelinesEnum, emptyInputs));
    assertEquals("pipelineInputs must be a JSON object", exception.getMessage());
  }

  @Test
  void validateInputsMultipleMissing() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object missingRequiredInputs =
        new LinkedHashMap<String, Object>(Map.of("not_an_input", "who cares"));
    ValidationException exception =
        assertThrows(
            ValidationException.class,
            () -> pipelinesService.validateInputs(pipelinesEnum, missingRequiredInputs));
    assertEquals(
        "Problems with pipelineInputs: multi_sample_vcf is required; new_integer_input is required",
        exception.getMessage());
  }

  @Test
  void validateInputsOneMissingOneWrongType() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a pipeline input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "Integer", true);

    pipelineInputDefinitionsRepository.save(newInput);

    Object missingRequiredInputs =
        new LinkedHashMap<String, Object>(Map.of("new_integer_input", "this is not an integer"));
    ValidationException exception =
        assertThrows(
            ValidationException.class,
            () -> pipelinesService.validateInputs(pipelinesEnum, missingRequiredInputs));
    assertEquals(
        "Problems with pipelineInputs: multi_sample_vcf is required; new_integer_input must be of type Integer",
        exception.getMessage());
  }

  @Test
  void validateInputsNotAMap() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Object inputs = new ArrayList<>(List.of("this", "is", "not", "a", "map"));
    assertThrows(
        ValidationException.class, () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
  }

  @Test
  void validateRequiredInputPresent() {
    PipelineInputDefinition requiredInput =
        new PipelineInputDefinition(1L, "required_input", "Integer", true);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(requiredInput));

    LinkedHashMap<String, Object> inputs =
        new LinkedHashMap<>(Map.of("required_input", "this is a string"));

    // even though it's the wrong type, the required input is present, so should not return an error
    // message
    assertTrue(pipelinesService.validateRequiredInputs(inputDefinitions, inputs).isEmpty());
  }

  private static Stream<Arguments> inputTypeValidationsValid() {
    return Stream.of(
        // arguments are: type specification, value to be tested, isRequired
        arguments("Integer", 123, true),
        arguments("Integer", "123", true),
        arguments("String", "I am a string", true),
        arguments("VCF", "path/to/file.vcf.gz", true),
        arguments("Array_String", List.of("this", "is", "a", "list", "of", "strings"), true),
        arguments(
            "Array_String",
            List.of(1, 2, 3),
            true), // including to remind us that this currently validates
        arguments("Array_VCF", List.of("path/to/file.vcf.gz"), true),
        arguments("String", "I am a string", false));
  }

  @ParameterizedTest
  @MethodSource("inputTypeValidationsValid")
  void validateInputTypeValid(String type, Object inputValue, Boolean isRequired) {
    PipelineInputDefinition requiredInput =
        new PipelineInputDefinition(1L, "input_name", type, isRequired);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(requiredInput));

    LinkedHashMap<String, Object> inputs = new LinkedHashMap<>(Map.of("input_name", inputValue));

    assertTrue(pipelinesService.validateInputTypes(inputDefinitions, inputs).isEmpty());
  }

  private static Stream<Arguments> inputTypeValidationsNotValid() {
    String commonTypeErrorMessage = "input_name must be of type %s";
    String vcfTypeErrorMessage = "input_name must be a path to a VCF file ending in .vcf.gz";
    String arrayStringTypeErrorMessage = "input_name must be an array of non-null strings";
    String arrayVCFTypeErrorMessage =
        "input_name must be an array of non-null paths to VCF files ending in .vcf.gz";
    String notNullErrorMessage = "input_name must not be null";

    String[] stringArrayWithNull = new String[] {null, "has a null value"};

    return Stream.of(
        // arguments are: type specification, value to be tested, isRequired, expected errorMessage
        // basic type checks
        arguments("Integer", "I am a string", true, commonTypeErrorMessage),
        arguments("Integer", 2.3, true, commonTypeErrorMessage),
        arguments("Integer", "2.3", true, commonTypeErrorMessage),
        arguments("String", "", true, "input_name must not be empty"),
        arguments(
            "String", List.of("this", "is", "not", "a", "string"), true, commonTypeErrorMessage),
        arguments("VCF", "path/to/file.vcf", true, vcfTypeErrorMessage),
        arguments("VCF", 3, true, vcfTypeErrorMessage),
        arguments(
            "Array_VCF",
            new String[] {"path/to/file.vcf.gz", "not a path"},
            true,
            arrayVCFTypeErrorMessage),
        arguments("Array_String", "I am a string", true, commonTypeErrorMessage),
        arguments(
            "Array_String",
            List.of("", "there's an empty string"),
            true,
            arrayStringTypeErrorMessage),
        arguments(
            "Array_String",
            new String[] {null, "list with null"},
            true,
            arrayStringTypeErrorMessage),
        arguments("Integer", "I am a string", false, commonTypeErrorMessage),
        // null checks
        arguments("String", null, false, notNullErrorMessage),
        arguments("Integer", null, true, notNullErrorMessage),
        arguments("Array_string", null, true, notNullErrorMessage),
        arguments("Array_String", stringArrayWithNull, true, arrayStringTypeErrorMessage),
        arguments("Array_VCF", null, true, notNullErrorMessage),
        arguments(
            "Array_VCF",
            new String[] {null, "list/with/null.vcf.gz"},
            true,
            arrayVCFTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("inputTypeValidationsNotValid")
  void validateInputTypeNotValid(
      String type, Object inputValue, Boolean isRequired, String errorMessage) {
    PipelineInputDefinition requiredInput =
        new PipelineInputDefinition(1L, "input_name", type, isRequired);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(requiredInput));

    LinkedHashMap<String, Object> inputs = new LinkedHashMap<>();
    inputs.put("input_name", inputValue);

    assertEquals(
        errorMessage.formatted(type),
        pipelinesService.validateInputTypes(inputDefinitions, inputs).get(0));
  }

  @Test
  void validateOptionalInputTypeValidAbsent() {
    PipelineInputDefinition requiredInput =
        new PipelineInputDefinition(1L, "required_input", "String", false);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(requiredInput));

    LinkedHashMap<String, Object> inputs = new LinkedHashMap<>();

    assertTrue(pipelinesService.validateInputTypes(inputDefinitions, inputs).isEmpty());
  }
}
