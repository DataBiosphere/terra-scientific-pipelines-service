package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.OutputDestination;
import bio.terra.cbas.model.ParameterDefinition;
import bio.terra.cbas.model.ParameterTypeDefinition;
import bio.terra.cbas.model.ParameterTypeDefinitionArray;
import bio.terra.cbas.model.ParameterTypeDefinitionPrimitive;
import bio.terra.cbas.model.PrimitiveParameterValueType;
import bio.terra.cbas.model.WorkflowInputDefinition;
import bio.terra.cbas.model.WorkflowOutputDefinition;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.ConstraintViolationException;
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
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class PipelinesServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelinesService pipelinesService;
  @Autowired PipelinesRepository pipelinesRepository;
  @MockBean SamService samService;
  @MockBean RawlsService rawlsService;

  @Test
  void getCorrectNumberOfPipelines() {
    // migrations insert one pipeline (imputation) so make sure we find it
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    UUID workspaceId = UUID.randomUUID();
    String workspaceBillingProject = "testTerraProject";
    String workspaceName = "testTerraWorkspaceName";
    String workspaceStorageContainerName = "testWorkspaceStorageContainerUrl";
    String workspaceGoogleProject = "testWorkspaceGoogleProject";

    // save a new version of the same pipeline
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.IMPUTATION_BEAGLE,
            1,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "wdlMethodName",
            "1.2.1",
            workspaceId,
            workspaceBillingProject,
            workspaceName,
            workspaceStorageContainerName,
            workspaceGoogleProject,
            null,
            null));

    pipelineList = pipelinesService.getPipelines();
    assertEquals(2, pipelineList.size());
    Pipeline savedPipeline = pipelineList.get(1);
    assertEquals(PipelinesEnum.IMPUTATION_BEAGLE, savedPipeline.getName());
    assertEquals(1, savedPipeline.getVersion());
    assertEquals("pipelineDisplayName", savedPipeline.getDisplayName());
    assertEquals("description", savedPipeline.getDescription());
    assertEquals("pipelineType", savedPipeline.getPipelineType());
    assertEquals("wdlUrl", savedPipeline.getWdlUrl());
    assertEquals("wdlMethodName", savedPipeline.getWdlMethodName());
    assertEquals("1.2.1", savedPipeline.getWdlMethodVersion());
    assertEquals(workspaceId, savedPipeline.getWorkspaceId());
    assertEquals(workspaceBillingProject, savedPipeline.getWorkspaceBillingProject());
    assertEquals(workspaceName, savedPipeline.getWorkspaceName());
    assertEquals(workspaceStorageContainerName, savedPipeline.getWorkspaceStorageContainerName());
    assertEquals(workspaceGoogleProject, savedPipeline.getWorkspaceGoogleProject());
  }

  @Test
  void getPipelineById() {
    PipelinesEnum imputationPipeline = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline p = pipelinesService.getPipeline(imputationPipeline);
    Long pipelineId = p.getId();

    Pipeline pById = pipelinesService.getPipelineById(pipelineId);
    assertEquals(p, pById);

    assertThrows(NotFoundException.class, () -> pipelinesService.getPipelineById(999L));
  }

  @Test
  void getPipelineInputDefinitions() {
    PipelinesEnum imputationPipeline = PipelinesEnum.IMPUTATION_BEAGLE;
    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipelinesService.getPipeline(imputationPipeline).getPipelineInputDefinitions();
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
              "Pipeline[pipelineName=%s, version=%s, displayName=%s, description=%s, pipelineType=%s, wdlUrl=%s, wdlMethodName=%s, wdlMethodVersion=%s, workspaceId=%s, workspaceBillingProject=%s, workspaceName=%s, workspaceStorageContainerName=%s, workspaceGoogleProject=%s]",
              p.getName(),
              p.getVersion(),
              p.getDisplayName(),
              p.getDescription(),
              p.getPipelineType(),
              p.getWdlUrl(),
              p.getWdlMethodName(),
              p.getWdlMethodVersion(),
              p.getWorkspaceId(),
              p.getWorkspaceBillingProject(),
              p.getWorkspaceName(),
              p.getWorkspaceStorageContainerName(),
              p.getWorkspaceGoogleProject()),
          p.toString());
    }
  }

  @Test
  void testPipelineHashCode() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      // 17 and 31 are hardcoded in this hashCode method of this class
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
              .append(p.getWdlMethodVersion())
              .append(p.getWorkspaceId())
              .append(p.getWorkspaceBillingProject())
              .append(p.getWorkspaceName())
              .append(p.getWorkspaceStorageContainerName())
              .append(p.getWorkspaceGoogleProject())
              .toHashCode(),
          p.hashCode());
    }
  }

  @Test
  void updatePipelineWorkspace() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    String newWdlMethodVersion = "0.13.1";

    String newWorkspaceStorageContainerName = "newTestWorkspaceStorageContainerUrl";
    String newWorkspaceGoogleProject = "newTestWorkspaceGoogleProject";
    WorkspaceDetails workspaceDetails =
        new WorkspaceDetails()
            .bucketName(newWorkspaceStorageContainerName)
            .googleProject(newWorkspaceGoogleProject);

    // make sure the current pipeline does not have the workspace info we're trying to update with
    assertNotEquals(newWorkspaceBillingProject, p.getWorkspaceBillingProject());
    assertNotEquals(newWorkspaceName, p.getWorkspaceName());
    assertNotEquals(newWorkspaceStorageContainerName, p.getWorkspaceStorageContainerName());
    assertNotEquals(newWorkspaceGoogleProject, p.getWorkspaceGoogleProject());
    assertNotEquals(newWdlMethodVersion, p.getWdlMethodVersion());

    // mocks
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("fakeToken");
    when(rawlsService.getWorkspaceDetails(
            "fakeToken", newWorkspaceBillingProject, newWorkspaceName))
        .thenReturn(workspaceDetails);
    when(rawlsService.getWorkspaceBucketName(workspaceDetails))
        .thenReturn(newWorkspaceStorageContainerName);
    when(rawlsService.getWorkspaceGoogleProject(workspaceDetails))
        .thenReturn(newWorkspaceGoogleProject);

    // update pipeline workspace id
    pipelinesService.updatePipelineWorkspace(
        pipelinesEnum, newWorkspaceBillingProject, newWorkspaceName, newWdlMethodVersion);
    p = pipelinesService.getPipeline(pipelinesEnum);

    // assert the workspace info has been updated
    assertEquals(newWorkspaceBillingProject, p.getWorkspaceBillingProject());
    assertEquals(newWorkspaceName, p.getWorkspaceName());
    assertEquals(newWdlMethodVersion, p.getWdlMethodVersion());

    // assert the fetched info from Rawls has been written
    assertEquals(newWorkspaceStorageContainerName, p.getWorkspaceStorageContainerName());
    assertEquals(newWorkspaceGoogleProject, p.getWorkspaceGoogleProject());

    // update pipeline with warp-like wdlMethodVersion
    String warpListWdlMethodVersion = "ImputationBeagle_development_v0.0.0";
    pipelinesService.updatePipelineWorkspace(
        pipelinesEnum, newWorkspaceBillingProject, newWorkspaceName, warpListWdlMethodVersion);

    p = pipelinesService.getPipeline(pipelinesEnum);

    assertEquals(warpListWdlMethodVersion, p.getWdlMethodVersion());
  }

  private static Stream<Arguments> badWdlMethodVersions() {
    return Stream.of(
        arguments("1.13.1"), // current imputation beagle pipeline is version 0 so this should fail
        arguments("blah.3.2"),
        arguments("0.2"),
        arguments("1.bhmm.2"),
        arguments("1.4.ok"),
        arguments("1.v3.4"),
        arguments(
            "ImputationBeagle_development_0.0.0"), // missing a "v" before the semantic version
        arguments("hiiv.1.4"));
  }

  @ParameterizedTest
  @MethodSource("badWdlMethodVersions")
  void updatePipelineWorkspaceBadWdlMethodVersion(String badWdlMethodVersion) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    assertThrows(
        ValidationException.class,
        () ->
            pipelinesService.updatePipelineWorkspace(
                pipelinesEnum, newWorkspaceBillingProject, newWorkspaceName, badWdlMethodVersion));
  }

  @Test
  void updatePipelineWorkspaceNullValueThrows() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum);

    String newWorkspaceName = "newTestTerraWorkspaceName";

    // make sure the current pipeline does not have the workspace info we're trying to update with
    assertNotEquals(newWorkspaceName, p.getWorkspaceName());

    // attempt to update pipeline info including null values should fail
    assertThrows(
        ConstraintViolationException.class,
        () ->
            pipelinesService.updatePipelineWorkspace(pipelinesEnum, null, newWorkspaceName, null));
  }

  static final String REQUIRED_STRING_INPUT_NAME = "outputBasename";
  static final String REQUIRED_VCF_INPUT_NAME = "multiSampleVcf";

  // input validation tests
  private static Stream<Arguments> inputValidations() {
    return Stream.of(
        // arguments: inputs, shouldPassValidation, list of strings expectedErrorMessage should
        // contain
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
            List.of(
                "Problem(s) with pipelineInputs:",
                "%s is required".formatted(REQUIRED_VCF_INPUT_NAME),
                "%s is required".formatted(REQUIRED_STRING_INPUT_NAME))),
        arguments(
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_STRING_INPUT_NAME, Arrays.asList("this is an array, not a string"))),
            false,
            List.of(
                "Problem(s) with pipelineInputs:",
                "%s is required".formatted(REQUIRED_VCF_INPUT_NAME),
                "%s must be a string".formatted(REQUIRED_STRING_INPUT_NAME))));
  }

  @ParameterizedTest
  @MethodSource("inputValidations")
  void validateInputs(
      Map<String, Object> inputs,
      Boolean shouldPassValidation,
      List<String> expectedErrorMessageStrings) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    List<PipelineInputDefinition> allInputDefinitions =
        pipelinesService.getPipeline(pipelinesEnum).getPipelineInputDefinitions();

    if (shouldPassValidation) {
      assertDoesNotThrow(
          () -> pipelinesService.validateUserProvidedInputs(allInputDefinitions, inputs));
    } else {
      ValidationException exception =
          assertThrows(
              ValidationException.class,
              () -> pipelinesService.validateUserProvidedInputs(allInputDefinitions, inputs));
      // this allows us to not care about the order in which the error messages are returned,
      // which depends on which pipelineInputDefinition was last updated in the db
      for (String expectedErrorMessage : expectedErrorMessageStrings) {
        assertTrue(exception.getMessage().contains(expectedErrorMessage));
      }
    }
  }

  private static Stream<Arguments> inputRequiredValidations() {
    return Stream.of(
        // arguments: isRequired, inputs, shouldPassValidation
        arguments(true, new HashMap<String, Object>(Map.of("inputName", "value")), true),
        arguments(true, new HashMap<String, Object>(), false),
        arguments(false, new HashMap<String, Object>(Map.of("inputName", "value")), true),
        arguments(false, new HashMap<String, Object>(), true));
  }

  @ParameterizedTest
  @MethodSource("inputRequiredValidations")
  void validateRequiredInputPresent(
      Boolean isRequired, Map<String, Object> inputs, Boolean shouldPassValidation) {
    PipelineInputDefinition inputDefinition =
        // note that pipelineId here is arbitrary since it's not used in the validate method
        new PipelineInputDefinition(
            1L,
            "inputName",
            "input_name",
            PipelineVariableTypesEnum.INTEGER,
            null,
            isRequired,
            true,
            null);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    if (shouldPassValidation) {
      assertTrue(pipelinesService.validateRequiredInputs(inputDefinitions, inputs).isEmpty());
    } else {
      assertEquals(1, pipelinesService.validateRequiredInputs(inputDefinitions, inputs).size());
      assertEquals(
          "inputName is required",
          pipelinesService.validateRequiredInputs(inputDefinitions, inputs).get(0));
    }
  }

  private static Stream<Arguments> inputTypeValidations() {
    return Stream.of(
        // arguments: type specification, file suffix if FILE type, value to be tested, whether it
        // should pass validation
        // INTEGER
        arguments(PipelineVariableTypesEnum.INTEGER, null, 123, true),
        arguments(PipelineVariableTypesEnum.INTEGER, null, "123", true),
        arguments(PipelineVariableTypesEnum.INTEGER, null, "I am a string", false),
        arguments(PipelineVariableTypesEnum.INTEGER, null, 2.3, false),
        arguments(PipelineVariableTypesEnum.INTEGER, null, "2.3", false),
        arguments(PipelineVariableTypesEnum.INTEGER, null, null, false),
        arguments(PipelineVariableTypesEnum.INTEGER, null, "", false),

        // STRING
        arguments(PipelineVariableTypesEnum.STRING, null, "I am a string", true),
        arguments(PipelineVariableTypesEnum.STRING, null, "  I am a string  ", true),
        arguments(
            PipelineVariableTypesEnum.STRING,
            null,
            List.of("this", "is", "not", "a", "string"),
            false),
        arguments(PipelineVariableTypesEnum.STRING, null, 123, false),
        arguments(PipelineVariableTypesEnum.STRING, null, null, false),
        arguments(PipelineVariableTypesEnum.STRING, null, "", false),

        // FILE
        arguments(PipelineVariableTypesEnum.FILE, ".vcf.gz", "path/to/file.vcf.gz", true),
        arguments(PipelineVariableTypesEnum.FILE, ".bed", "path/to/file.bed", true),
        arguments(PipelineVariableTypesEnum.FILE, ".vcf.gz", "path/to/file.vcf", false),
        arguments(PipelineVariableTypesEnum.FILE, ".vcf.gz", 3, false),
        arguments(PipelineVariableTypesEnum.FILE, ".vcf.gz", null, false),
        arguments(PipelineVariableTypesEnum.FILE, ".vcf.gz", "", false),

        // STRING_ARRAY
        arguments(
            PipelineVariableTypesEnum.STRING_ARRAY,
            null,
            Arrays.asList("this", "is", "a", "list", "of", "strings"),
            true),
        arguments(
            PipelineVariableTypesEnum.STRING_ARRAY,
            null,
            Arrays.asList("this ", " is", " a ", "list  ", "  of", "  strings  "),
            true),
        arguments(PipelineVariableTypesEnum.STRING_ARRAY, null, "I am not an array", false),
        arguments(PipelineVariableTypesEnum.STRING_ARRAY, null, Arrays.asList(1, 2, 3), false),
        arguments(PipelineVariableTypesEnum.STRING_ARRAY, null, null, false),
        arguments(PipelineVariableTypesEnum.STRING_ARRAY, null, List.of(), false),
        arguments(
            PipelineVariableTypesEnum.STRING_ARRAY,
            null,
            Arrays.asList("", "array with empty string"),
            false),
        arguments(
            PipelineVariableTypesEnum.STRING_ARRAY,
            null,
            Arrays.asList(null, "array with null"),
            false),

        // FILE_ARRAY
        arguments(
            PipelineVariableTypesEnum.FILE_ARRAY, ".vcf.gz", List.of("path/to/file.vcf.gz"), true),
        arguments(
            PipelineVariableTypesEnum.FILE_ARRAY,
            ".vcf.gz",
            List.of(" path/to/file.vcf.gz  "),
            true),
        arguments(
            PipelineVariableTypesEnum.FILE_ARRAY, ".vcf.gz", "this/is/not/an/array.vcf.gz", false),
        arguments(
            PipelineVariableTypesEnum.FILE_ARRAY,
            ".vcf.gz",
            Arrays.asList("path/to/file.vcf.gz", "not a path"),
            false),
        arguments(PipelineVariableTypesEnum.FILE_ARRAY, ".vcf.gz", null, false),
        arguments(PipelineVariableTypesEnum.FILE_ARRAY, ".vcf.gz", List.of(), false),
        arguments(
            PipelineVariableTypesEnum.FILE_ARRAY,
            ".vcf.gz",
            Arrays.asList(null, "list/with/null.vcf.gz"),
            false));
  }

  @ParameterizedTest
  @MethodSource("inputTypeValidations")
  void validateInputType(
      PipelineVariableTypesEnum inputType,
      String fileSuffix,
      Object inputValue,
      Boolean shouldPassValidation) {
    PipelineInputDefinition inputDefinition =
        new PipelineInputDefinition(
            1L, "inputName", "input_name", inputType, fileSuffix, true, true, null);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    Map<String, Object> inputs = new HashMap<>();
    inputs.put("inputName", inputValue);

    // error message contents are tested in PipelineInputTypesEnumTest
    assertEquals(
        shouldPassValidation,
        pipelinesService.validateInputTypes(inputDefinitions, inputs).isEmpty());
  }

  @Test
  void constructImputationInputsSuccess() {
    Map<String, Object> userProvidedInputs = TestUtils.TEST_PIPELINE_INPUTS;

    PipelinesEnum pipelineEnum = PipelinesEnum.IMPUTATION_BEAGLE;
    Pipeline pipeline = pipelinesRepository.findByName(pipelineEnum);
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

  @Test
  void extractUserProvidedFileInputNames() {
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>();
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input1", "input_1", PipelineVariableTypesEnum.STRING, null, true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input2", "input_2", PipelineVariableTypesEnum.INTEGER, null, false, true, "1"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input3",
            "input_3",
            PipelineVariableTypesEnum.FILE,
            ".vcf.gz",
            true,
            false,
            "not/a/user/provided/input.vcf.gz"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input4", "input_4", PipelineVariableTypesEnum.FILE, ".vcf.gz", true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input5",
            "input_5",
            PipelineVariableTypesEnum.FILE_ARRAY,
            ".vcf.gz",
            true,
            true,
            null));

    List<String> userProvidedFileInputNames =
        pipelinesService.extractUserProvidedFileInputNames(inputDefinitions);

    // the only user-provided inputs that we currently recognize as files are VCFs, i.e. input4
    assertEquals(1, userProvidedFileInputNames.size());
    assertTrue(userProvidedFileInputNames.contains("input4"));
  }

  private static Stream<Arguments> mapVariableTypeToCbasParameterTypeArguments() {
    ParameterTypeDefinition stringParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.STRING)
            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
    ParameterTypeDefinition fileParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.FILE)
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
    ParameterTypeDefinition fileArrayParameterTypeResponse =
        new ParameterTypeDefinitionArray()
            .nonEmpty(true)
            .arrayType(
                new ParameterTypeDefinitionPrimitive()
                    .primitiveType(PrimitiveParameterValueType.FILE)
                    .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
            .type(ParameterTypeDefinition.TypeEnum.ARRAY);
    return Stream.of(
        // arguments: type specification, expected response
        arguments(PipelineVariableTypesEnum.STRING, stringParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.FILE, fileParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.INTEGER, integerParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.STRING_ARRAY, stringArrayParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.FILE_ARRAY, fileArrayParameterTypeResponse));
  }

  @ParameterizedTest
  @MethodSource("mapVariableTypeToCbasParameterTypeArguments")
  void mapVariableTypeToCbasParameterType(
      PipelineVariableTypesEnum inputType, ParameterTypeDefinition expectedResponse) {
    assertEquals(expectedResponse, pipelinesService.mapVariableTypeToCbasParameterType(inputType));
  }

  @Test
  void prepareCbasWorkflowInputRecordLookupDefinitions() {
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>();
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input1", "input_1", PipelineVariableTypesEnum.STRING, null, true, true, null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L, "input2", "input_2", PipelineVariableTypesEnum.INTEGER, null, false, true, "1"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input3",
            "input_3",
            PipelineVariableTypesEnum.STRING_ARRAY,
            null,
            true,
            false,
            "[\"1\", \"2\"]"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input4",
            "input_4",
            PipelineVariableTypesEnum.FILE,
            ".vcf.gz",
            false,
            false,
            "fake/file.vcf.gz"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input5",
            "input_5",
            PipelineVariableTypesEnum.FILE_ARRAY,
            ".vcf.gz",
            true,
            true,
            null));

    String testWdlName = "aFakeWdl";

    List<WorkflowInputDefinition> cbasWorkflowInputDefinitions =
        pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(
            inputDefinitions, testWdlName);

    assertEquals(inputDefinitions.size(), cbasWorkflowInputDefinitions.size());
    for (int i = 0; i < inputDefinitions.size(); i++) {
      // the input type should be the object returned by the mapVariableTypeToCbasParameterType
      // method
      assertEquals(
          pipelinesService.mapVariableTypeToCbasParameterType(inputDefinitions.get(i).getType()),
          cbasWorkflowInputDefinitions.get(i).getInputType());
      // the input name should be the wdl name concatenated with the input name
      assertEquals(
          "%s.%s".formatted(testWdlName, inputDefinitions.get(i).getWdlVariableName()),
          cbasWorkflowInputDefinitions.get(i).getInputName());
      // the source should be a record lookup
      assertEquals(
          ParameterDefinition.TypeEnum.RECORD_LOOKUP,
          cbasWorkflowInputDefinitions.get(i).getSource().getType());
    }
  }

  @Test
  void prepareCbasWorkflowOutputRecordUpdateDefinitions() {
    List<PipelineOutputDefinition> outputDefinitions = new ArrayList<>();
    outputDefinitions.add(
        new PipelineOutputDefinition(1L, "output1", "output_1", PipelineVariableTypesEnum.FILE));
    outputDefinitions.add(
        new PipelineOutputDefinition(1L, "output2", "output_2", PipelineVariableTypesEnum.STRING));
    String testWdlName = "aFakeWdl";

    List<WorkflowOutputDefinition> cbasWorkflowOutputDefinitions =
        pipelinesService.prepareCbasWorkflowOutputRecordUpdateDefinitions(
            outputDefinitions, testWdlName);

    assertEquals(outputDefinitions.size(), cbasWorkflowOutputDefinitions.size());
    for (int i = 0; i < outputDefinitions.size(); i++) {
      // the output name should be the wdl name concatenated with the output name
      assertEquals(
          "%s.%s".formatted(testWdlName, outputDefinitions.get(i).getWdlVariableName()),
          cbasWorkflowOutputDefinitions.get(i).getOutputName());
      // the output type should be the object returned by the mapVariableTypeToCbasParameterType
      // method
      assertEquals(
          pipelinesService.mapVariableTypeToCbasParameterType(outputDefinitions.get(i).getType()),
          cbasWorkflowOutputDefinitions.get(i).getOutputType());
      // the destination type should be a record update
      assertEquals(
          OutputDestination.TypeEnum.RECORD_UPDATE,
          cbasWorkflowOutputDefinitions.get(i).getDestination().getType());
    }
  }
}
