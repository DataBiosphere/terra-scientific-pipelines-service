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
      assertDoesNotThrow(() -> PipelineInputTypesEnum.valueOf(p.getType()));
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
    newInput.setType("INTEGER");
    newInput.setIsRequired(false);

    pipelineInputDefinitionsRepository.save(newInput);

    pipeline = pipelinesRepository.findByName(PipelinesEnum.IMPUTATION_BEAGLE.getValue());
    pipelineInputDefinitions = pipeline.getPipelineInputDefinitions();
    assertEquals(2, pipelineInputDefinitions.size());

    PipelineInputDefinition savedInput = pipelineInputDefinitions.get(1);
    assertEquals("newInput", savedInput.getName());
    assertEquals("INTEGER", savedInput.getType());
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
  private static Stream<Arguments> inputValidations() {
    return Stream.of(
        // arguments: isRequired, inputs, shouldPassValidation, expectedErrorMessage
        arguments(
            new LinkedHashMap<String, Object>(
                Map.of("multi_sample_vcf", "this/is/a/vcf/path.vcf.gz", "new_integer_input", 123)),
            true,
            null),
        arguments(new Object(), false, "pipelineInputs must be a JSON object"),
        arguments(
            new ArrayList<>(List.of("this", "is", "not", "a", "map")),
            false,
            "pipelineInputs must be a JSON object"),
        arguments(
            new LinkedHashMap<String, Object>(Map.of("not_an_input", "who cares")),
            false,
            "Problem(s) with pipelineInputs: multi_sample_vcf is required; new_integer_input is required"),
        arguments(
            new LinkedHashMap<String, Object>(
                Map.of("new_integer_input", "this is not an integer")),
            false,
            "Problem(s) with pipelineInputs: multi_sample_vcf is required; new_integer_input must be of type INTEGER"));
  }

  @ParameterizedTest
  @MethodSource("inputValidations")
  void validateInputs(Object inputs, Boolean shouldPassValidation, String expectedErrorMessage) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelinesEnum.getValue());

    // add a required INTEGER input to the imputation pipeline
    PipelineInputDefinition newInput =
        new PipelineInputDefinition(pipeline.getId(), "new_integer_input", "INTEGER", true);

    pipelineInputDefinitionsRepository.save(newInput);

    if (shouldPassValidation) {
      assertDoesNotThrow(() -> pipelinesService.validateInputs(pipelinesEnum, inputs));
    } else {
      ValidationException exception =
          assertThrows(
              ValidationException.class,
              () -> pipelinesService.validateInputs(pipelinesEnum, inputs));
      assertEquals(expectedErrorMessage, exception.getMessage());
    }
  }

  private static Stream<Arguments> inputRequiredValidations() {
    return Stream.of(
        // arguments: isRequired, inputs, shouldPassValidation
        arguments(true, new LinkedHashMap<String, Object>(Map.of("input_name", "value")), true),
        arguments(true, new LinkedHashMap<String, Object>(), false),
        arguments(false, new LinkedHashMap<String, Object>(Map.of("input_name", "value")), true),
        arguments(false, new LinkedHashMap<String, Object>(), true));
  }

  @ParameterizedTest
  @MethodSource("inputRequiredValidations")
  void validateRequiredInputPresent(
      Boolean isRequired, Map<String, Object> inputs, Boolean shouldPassValidation) {
    PipelineInputDefinition inputDefinition =
        new PipelineInputDefinition(1L, "input_name", "INTEGER", isRequired);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    if (shouldPassValidation) {
      assertTrue(pipelinesService.validateRequiredInputs(inputDefinitions, inputs).isEmpty());
    } else {
      assertEquals(1, pipelinesService.validateRequiredInputs(inputDefinitions, inputs).size());
      assertEquals(
          "input_name is required",
          pipelinesService.validateRequiredInputs(inputDefinitions, inputs).get(0));
    }
  }

  private static Stream<Arguments> inputTypeValidations() {
    // error messages
    String commonTypeErrorMessage = "input_name must be of type %s";
    String vcfTypeErrorMessage = "input_name must be a path to a VCF file ending in .vcf.gz";
    String arrayStringTypeErrorMessage = "input_name must be an array of non-null strings";
    String arrayVCFTypeErrorMessage =
        "input_name must be an array of non-null paths to VCF files ending in .vcf.gz";
    String notNullErrorMessage = "input_name must not be null";
    String notEmptyErrorMessage = "input_name must not be empty";
    String emptyArrayErrorMessage = "input_name must not be an empty array";

    return Stream.of(
        // arguments: type specification, value to be tested, whether it should pass validation,
        // expected error message if not
        arguments("INTEGER", 123, true, null),
        arguments("INTEGER", "123", true, null),
        arguments("STRING", "I am a string", true, null),
        arguments(
            "STRING", 123, true, null), // including to remind us that this currently validates
        arguments("VCF", "path/to/file.vcf.gz", true, null),
        arguments("ARRAY_STRING", List.of("this", "is", "a", "list", "of", "strings"), true, null),
        arguments(
            "ARRAY_STRING",
            List.of(1, 2, 3),
            true,
            null), // including to remind us that this currently validates
        arguments("ARRAY_VCF", new String[] {"path/to/file.vcf.gz"}, true, null),
        // basic type checks that should fail validation (produce an error message)
        arguments(
            "STRING", List.of("this", "is", "not", "a", "string"), false, commonTypeErrorMessage),
        arguments("INTEGER", "I am a string", false, commonTypeErrorMessage),
        arguments("INTEGER", 2.3, false, commonTypeErrorMessage),
        arguments("INTEGER", "2.3", false, commonTypeErrorMessage),
        arguments("VCF", "path/to/file.vcf", false, vcfTypeErrorMessage),
        arguments("VCF", 3, false, vcfTypeErrorMessage),
        arguments("ARRAY_STRING", "I am not an array", false, commonTypeErrorMessage),
        arguments("ARRAY_VCF", "this/is/not/an/array.vcf.gz", false, commonTypeErrorMessage),
        arguments(
            "ARRAY_VCF",
            new String[] {"path/to/file.vcf.gz", "not a path"},
            false,
            arrayVCFTypeErrorMessage),
        // null and empty checks
        arguments("STRING", null, false, notNullErrorMessage),
        arguments("STRING", "", false, notEmptyErrorMessage),
        arguments("INTEGER", null, false, notNullErrorMessage),
        arguments("INTEGER", "", false, notEmptyErrorMessage),
        arguments("ARRAY_STRING", null, false, notNullErrorMessage),
        arguments("ARRAY_STRING", new String[] {}, false, emptyArrayErrorMessage),
        arguments(
            "ARRAY_STRING",
            List.of("", "array with empty string"),
            false,
            arrayStringTypeErrorMessage),
        arguments(
            "ARRAY_STRING",
            new String[] {null, "array with null"},
            false,
            arrayStringTypeErrorMessage),
        arguments("ARRAY_VCF", null, false, notNullErrorMessage),
        arguments("ARRAY_VCF", new String[] {}, false, emptyArrayErrorMessage),
        arguments(
            "ARRAY_VCF",
            new String[] {null, "list/with/null.vcf.gz"},
            false,
            arrayVCFTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("inputTypeValidations")
  void validateInputType(
      String inputType,
      Object inputValue,
      Boolean shouldPassValidation,
      String expectedErrorMessage) {
    PipelineInputDefinition inputDefinition =
        new PipelineInputDefinition(1L, "input_name", inputType, true);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    LinkedHashMap<String, Object> inputs = new LinkedHashMap<>();
    inputs.put("input_name", inputValue);

    if (shouldPassValidation) {
      assertTrue(pipelinesService.validateInputTypes(inputDefinitions, inputs).isEmpty());
    } else {
      assertEquals(1, pipelinesService.validateInputTypes(inputDefinitions, inputs).size());
      assertEquals(
          String.format(expectedErrorMessage, inputType),
          pipelinesService.validateInputTypes(inputDefinitions, inputs).get(0));
    }
  }
}
