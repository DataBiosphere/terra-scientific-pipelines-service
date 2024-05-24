package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.cbas.model.ParameterDefinition;
import bio.terra.cbas.model.ParameterTypeDefinition;
import bio.terra.cbas.model.ParameterTypeDefinitionArray;
import bio.terra.cbas.model.ParameterTypeDefinitionPrimitive;
import bio.terra.cbas.model.PrimitiveParameterValueType;
import bio.terra.cbas.model.WorkflowInputDefinition;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.repositories.PipelineInputDefinitionsRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
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
  void getPipelineInputDefinitions() {
    PipelinesEnum imputationPipeline = PipelinesEnum.IMPUTATION_BEAGLE;
    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipelinesService.getPipeline(imputationPipeline).getPipelineInputDefinitions().stream()
            .collect(Collectors.toList());
    List<PipelineInputDefinition> userProvidedPipelineInputDefinitions =
        pipelinesService.extractUserProvidedInputDefinitions(allPipelineInputDefinitions);
    List<PipelineInputDefinition> serviceProvidedPipelineInputDefinitions =
        pipelinesService.extractServiceProvidedInputDefinitions(allPipelineInputDefinitions);

    assertEquals(
        allPipelineInputDefinitions.size(),
        userProvidedPipelineInputDefinitions.size()
            + serviceProvidedPipelineInputDefinitions.size());
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

  static final String REQUIRED_STRING_INPUT_NAME = "output_basename";
  static final String REQUIRED_VCF_INPUT_NAME = "multi_sample_vcf";

  // input validation tests
  private static Stream<Arguments> inputValidations() {
    return Stream.of(
        // arguments: inputs, shouldPassValidation, expectedErrorMessage
        arguments(
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_STRING_INPUT_NAME,
                    "the-basename-value-for-my-output",
                    REQUIRED_VCF_INPUT_NAME,
                    "this/is/a/vcf/path.vcf.gz")),
            true,
            null),
        arguments(
            new HashMap<String, Object>(Map.of("not_an_input", "who cares")),
            false,
            "Problem(s) with pipelineInputs: %s is required; %s is required"
                .formatted(REQUIRED_VCF_INPUT_NAME, REQUIRED_STRING_INPUT_NAME)),
        arguments(
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_STRING_INPUT_NAME, Arrays.asList("this is an array, not a string"))),
            false,
            "Problem(s) with pipelineInputs: %s is required; %s must be a string"
                .formatted(REQUIRED_VCF_INPUT_NAME, REQUIRED_STRING_INPUT_NAME)));
  }

  @ParameterizedTest
  @MethodSource("inputValidations")
  void validateInputs(
      Map<String, Object> inputs, Boolean shouldPassValidation, String expectedErrorMessage) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    List<PipelineInputDefinition> allInputDefinitions =
        pipelinesService.getPipeline(pipelinesEnum).getPipelineInputDefinitions().stream()
            .collect(Collectors.toList());

    if (shouldPassValidation) {
      assertDoesNotThrow(
          () -> pipelinesService.validateUserProvidedInputs(allInputDefinitions, inputs));
    } else {
      ValidationException exception =
          assertThrows(
              ValidationException.class,
              () -> pipelinesService.validateUserProvidedInputs(allInputDefinitions, inputs));
      assertEquals(expectedErrorMessage, exception.getMessage());
    }
  }

  private static Stream<Arguments> inputRequiredValidations() {
    return Stream.of(
        // arguments: isRequired, inputs, shouldPassValidation
        arguments(true, new HashMap<String, Object>(Map.of("input_name", "value")), true),
        arguments(true, new HashMap<String, Object>(), false),
        arguments(false, new HashMap<String, Object>(Map.of("input_name", "value")), true),
        arguments(false, new HashMap<String, Object>(), true));
  }

  @ParameterizedTest
  @MethodSource("inputRequiredValidations")
  void validateRequiredInputPresent(
      Boolean isRequired, Map<String, Object> inputs, Boolean shouldPassValidation) {
    PipelineInputDefinition inputDefinition =
        // note that pipelineId here is arbitrary since it's not used in the validate method
        new PipelineInputDefinition(
            1L, "input_name", PipelineInputTypesEnum.INTEGER, isRequired, true, null);
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
    return Stream.of(
        // arguments: type specification, value to be tested, whether it should pass validation
        // INTEGER
        arguments(PipelineInputTypesEnum.INTEGER, 123, true),
        arguments(PipelineInputTypesEnum.INTEGER, "123", true),
        arguments(PipelineInputTypesEnum.INTEGER, "I am a string", false),
        arguments(PipelineInputTypesEnum.INTEGER, 2.3, false),
        arguments(PipelineInputTypesEnum.INTEGER, "2.3", false),
        arguments(PipelineInputTypesEnum.INTEGER, null, false),
        arguments(PipelineInputTypesEnum.INTEGER, "", false),

        // STRING
        arguments(PipelineInputTypesEnum.STRING, "I am a string", true),
        arguments(PipelineInputTypesEnum.STRING, "  I am a string  ", true),
        arguments(
            PipelineInputTypesEnum.STRING, List.of("this", "is", "not", "a", "string"), false),
        arguments(PipelineInputTypesEnum.STRING, 123, false),
        arguments(PipelineInputTypesEnum.STRING, null, false),
        arguments(PipelineInputTypesEnum.STRING, "", false),

        // VCF
        arguments(PipelineInputTypesEnum.VCF, "path/to/file.vcf.gz", true),
        arguments(PipelineInputTypesEnum.VCF, "path/to/file.vcf", false),
        arguments(PipelineInputTypesEnum.VCF, 3, false),
        arguments(PipelineInputTypesEnum.VCF, null, false),
        arguments(PipelineInputTypesEnum.VCF, "", false),

        // STRING_ARRAY
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            Arrays.asList("this", "is", "a", "list", "of", "strings"),
            true),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            Arrays.asList("this ", " is", " a ", "list  ", "  of", "  strings  "),
            true),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, "I am not an array", false),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, Arrays.asList(1, 2, 3), false),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, null, false),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, List.of(), false),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            Arrays.asList("", "array with empty string"),
            false),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY, Arrays.asList(null, "array with null"), false),

        // VCF_ARRAY
        arguments(PipelineInputTypesEnum.VCF_ARRAY, List.of("path/to/file.vcf.gz"), true),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, List.of(" path/to/file.vcf.gz  "), true),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, "this/is/not/an/array.vcf.gz", false),
        arguments(
            PipelineInputTypesEnum.VCF_ARRAY,
            Arrays.asList("path/to/file.vcf.gz", "not a path"),
            false),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, null, false),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, List.of(), false),
        arguments(
            PipelineInputTypesEnum.VCF_ARRAY, Arrays.asList(null, "list/with/null.vcf.gz"), false));
  }

  @ParameterizedTest
  @MethodSource("inputTypeValidations")
  void validateInputType(
      PipelineInputTypesEnum inputType, Object inputValue, Boolean shouldPassValidation) {
    PipelineInputDefinition inputDefinition =
        new PipelineInputDefinition(1L, "input_name", inputType, true, true, null);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    Map<String, Object> inputs = new HashMap<>();
    inputs.put("input_name", inputValue);

    // error message contents are tested in PipelineInputTypesEnumTest
    assertEquals(
        shouldPassValidation,
        pipelinesService.validateInputTypes(inputDefinitions, inputs).isEmpty());
  }

  @Test
  void constructImputationInputsSuccess() {
    Map<String, Object> userProvidedInputs = TestUtils.TEST_PIPELINE_INPUTS;

    PipelinesEnum pipelineEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelineEnum.getValue());
    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipeline.getPipelineInputDefinitions();

    List<PipelineInputDefinition> serviceProvidedPipelineInputDefinitions =
        pipelinesService.extractServiceProvidedInputDefinitions(allPipelineInputDefinitions);

    // this should add the service-provided inputs to the one user-provided input in
    // testPipelineInputs
    Map<String, Object> allPipelineInputs =
        pipelinesService.constructRawInputs(allPipelineInputDefinitions, userProvidedInputs);

    Integer totalInputs =
        userProvidedInputs.size() + serviceProvidedPipelineInputDefinitions.size();

    assertNotNull(allPipelineInputs);
    for (String inputName : userProvidedInputs.keySet()) {
      assertTrue(allPipelineInputs.containsKey(inputName));
    }
    for (String inputName :
        serviceProvidedPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())) {
      assertTrue(allPipelineInputs.containsKey(inputName));
    }
    assertEquals(totalInputs, allPipelineInputs.size());
  }

  private static Stream<Arguments> mapInputTypeToCbasParameterTypeArguments() {
    ParameterTypeDefinition stringParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.STRING)
            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
    ParameterTypeDefinition integerParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.INT)
            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
    ParameterTypeDefinition stringArrayParameterTypeResponse =
        new ParameterTypeDefinitionArray()
            .nonEmpty(true)
            .arrayType(
                new ParameterTypeDefinitionPrimitive()
                    .primitiveType(PrimitiveParameterValueType.STRING)
                    .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
            .type(ParameterTypeDefinition.TypeEnum.ARRAY);
    return Stream.of(
        // arguments: type specification, expected response
        arguments(PipelineInputTypesEnum.STRING, stringParameterTypeResponse),
        arguments(PipelineInputTypesEnum.VCF, stringParameterTypeResponse),
        arguments(PipelineInputTypesEnum.INTEGER, integerParameterTypeResponse),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, stringArrayParameterTypeResponse),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, stringArrayParameterTypeResponse));
  }

  @ParameterizedTest
  @MethodSource("mapInputTypeToCbasParameterTypeArguments")
  void mapInputTypeToCbasParameterType(
      PipelineInputTypesEnum inputType, ParameterTypeDefinition expectedResponse) {
    assertEquals(expectedResponse, pipelinesService.mapInputTypeToCbasParameterType(inputType));
  }

  @Test
  void prepareCbasWorkflowInputRecordLookupDefinitions() {
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>();
    inputDefinitions.add(
        new PipelineInputDefinition(1L, "input1", PipelineInputTypesEnum.STRING, true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input2", PipelineInputTypesEnum.INTEGER, true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input3", PipelineInputTypesEnum.STRING_ARRAY, true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(1L, "input4", PipelineInputTypesEnum.VCF, true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input5", PipelineInputTypesEnum.VCF_ARRAY, true, true, null));

    String testWdlName = "aFakeWdl";

    List<WorkflowInputDefinition> cbasWorkflowInputDefinitions =
        pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(
            inputDefinitions, testWdlName);

    assertEquals(inputDefinitions.size(), cbasWorkflowInputDefinitions.size());
    for (int i = 0; i < inputDefinitions.size(); i++) {
      // the input type should be the object returned by the mapInputTypeToCbasParameterType method
      assertEquals(
          pipelinesService.mapInputTypeToCbasParameterType(inputDefinitions.get(i).getType()),
          cbasWorkflowInputDefinitions.get(i).getInputType());
      // the input name should be the wdl name concatenated with the input name
      assertEquals(
          "%s.%s".formatted(testWdlName, inputDefinitions.get(i).getName()),
          cbasWorkflowInputDefinitions.get(i).getInputName());
      // the source should be a record lookup
      assertEquals(
          ParameterDefinition.TypeEnum.RECORD_LOOKUP,
          cbasWorkflowInputDefinitions.get(i).getSource().getType());
    }
  }
}
