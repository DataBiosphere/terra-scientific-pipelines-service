package bio.terra.pipelines.service;

import static bio.terra.pipelines.testutils.TestUtils.createTestPipelineWithId;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.ArgumentMatchers.anyString;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.DownloadOutputCall;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.DownloadOutputCallRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutputSignedUrls;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.rawls.model.Entity;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class PipelineInputsOutputsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelineInputsOutputsService pipelineInputsOutputsService;

  @Autowired PipelinesService pipelinesService;
  @Autowired PipelineInputsRepository pipelineInputsRepository;
  @Autowired PipelineOutputsRepository pipelineOutputsRepository;
  @Autowired DownloadOutputCallRepository downloadOutputCallRepository;

  @Autowired PipelineRunsRepository pipelineRunsRepository;

  @MockitoBean private GcsService mockGcsService;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @Test
  void prepareFileInputs() throws MalformedURLException {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generatePutObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceGoogleProject()),
            eq(testPipelineWithId.getWorkspaceStorageContainerName()),
            anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineInputsOutputsService.prepareFileInputs(
            testPipelineWithId, testJobId, userPipelineInputs, false);

    assertEquals(userPipelineInputs.size(), formattedPipelineFileInputs.size());
    assertEquals(
        fakeUrl.toString(), formattedPipelineFileInputs.get(fileInputKeyName).get("signedUrl"));
    assertEquals(
        "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s' | cat"
            .formatted(fileInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(fileInputKeyName).get("curlCommand"));
  }

  @Test
  void prepareFileInputsResumable() throws MalformedURLException {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generateResumablePostObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceGoogleProject()),
            eq(testPipelineWithId.getWorkspaceStorageContainerName()),
            anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineInputsOutputsService.prepareFileInputs(
            testPipelineWithId, testJobId, userPipelineInputs, true);

    assertEquals(userPipelineInputs.size(), formattedPipelineFileInputs.size());
    assertEquals(
        fakeUrl.toString(), formattedPipelineFileInputs.get(fileInputKeyName).get("signedUrl"));
    assertEquals(
        "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s $(curl -s -i -X POST -H 'x-goog-resumable: start' '%s' | grep -i '^Location:' | cut -d' ' -f2- | tr -d '\r') | cat"
            .formatted(fileInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(fileInputKeyName).get("curlCommand"));
  }

  @Test
  void extractPipelineOutputsFromEntity() {
    // test that the method correctly extracts the outputs from the entity
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(
        Map.of("output_string", "string", "output_integer", 123, "output_boolean_optional", false));

    Map<String, Object> extractedOutputs =
        pipelineInputsOutputsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);

    assertEquals(3, extractedOutputs.size());
    // the method should also have converted the wdlVariableName key to the camelCase outputName key
    assertEquals("string", extractedOutputs.get("outputString"));
    assertEquals(123, extractedOutputs.get("outputInteger"));
    assertEquals(false, extractedOutputs.get("outputBooleanOptional"));
  }

  @Test
  void extractPipelineOutputsFromEntityMissingOutput() {
    // test that the method correctly throws an error if a required output is missing
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(
        Map.of("outputInteger", 123, "outputBoolean", false)); // missing outputString

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineInputsOutputsService.extractPipelineOutputsFromEntity(
                outputDefinitions, entity));
  }

  @Test
  void extractPipelineOutputsFromEntityEmptyRequiredOutput() {
    // test that the method correctly throws an error if a required output is empty
    List<PipelineOutputDefinition> outputDefinitions =
        TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
    Entity entity = new Entity();
    entity.setAttributes(Map.of("outputString", "", "outputInteger", 123, "outputBoolean", false));

    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineInputsOutputsService.extractPipelineOutputsFromEntity(
                outputDefinitions, entity));
  }

  @Test
  void extractPipelineOutputsFromEntityEmptyOptionalOutput() {
    // test that the method succeeds if an optional output is empty
    List<PipelineOutputDefinition> outputDefinitions =
        new ArrayList<>(
            List.of(
                new PipelineOutputDefinition(
                    3L,
                    "outputString",
                    "output_string",
                    "output string",
                    "description",
                    PipelineVariableTypesEnum.STRING,
                    false)));
    Entity entity = new Entity();
    entity.setAttributes(Map.of("outputString", ""));

    Map<String, Object> extractedOutputs =
        pipelineInputsOutputsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);

    assertEquals(1, extractedOutputs.size());
    // the method should also have converted the wdlVariableName key to the camelCase outputName key
    assertNull(extractedOutputs.get("outputString"));
  }

  @Test
  void extractPipelineOutputsFromEntityMissingOptionalOutput() {
    // test that the method succeeds if an optional output is missing
    List<PipelineOutputDefinition> outputDefinitions =
        new ArrayList<>(
            List.of(
                new PipelineOutputDefinition(
                    3L,
                    "outputString",
                    "output_string",
                    "output string",
                    "description",
                    PipelineVariableTypesEnum.STRING,
                    false)));
    Entity entity = new Entity();
    entity.setAttributes(Map.of()); // missing outputString

    Map<String, Object> extractedOutputs =
        pipelineInputsOutputsService.extractPipelineOutputsFromEntity(outputDefinitions, entity);

    assertEquals(1, extractedOutputs.size());
    // the method should also have converted the wdlVariableName key to the camelCase outputName key
    assertNull(extractedOutputs.get("outputString"));
  }

  @Test
  void getPipelineRunOutputs() {
    // define pipeline with outputs
    Pipeline testPipeline = createTestPipelineWithId();
    testPipeline.setPipelineOutputDefinitions(TestUtils.TEST_PIPELINE_OUTPUT_DEFINITIONS_WITH_FILE);

    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(testJobId);
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRun.setPipeline(testPipeline);
    pipelineRunsRepository.save(pipelineRun);

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(pipelineRun.getId());
    pipelineOutput.setOutputs(
        pipelineInputsOutputsService.mapToString(TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE));
    pipelineOutputsRepository.save(pipelineOutput);

    Map<String, Object> retrievedOutputs =
        pipelineInputsOutputsService.getPipelineRunOutputs(pipelineRun);

    assertEquals(TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.size(), retrievedOutputs.size());
    for (String outputKey : TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.keySet()) {
      assertEquals(
          TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED.get(outputKey),
          retrievedOutputs.get(outputKey));
    }
  }

  @Test
  void generatePipelineRunOutputSignedUrls() throws MalformedURLException {
    // define pipeline with outputs: 1 file output and 1 string output
    Pipeline testPipeline = createTestPipelineWithId();
    testPipeline.setPipelineOutputDefinitions(TestUtils.TEST_PIPELINE_OUTPUT_DEFINITIONS_WITH_FILE);

    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(testJobId);
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRun.setPipeline(testPipeline);
    pipelineRunsRepository.save(pipelineRun);

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(pipelineRun.getId());
    pipelineOutput.setOutputs(
        pipelineInputsOutputsService.mapToString(TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE));
    pipelineOutputsRepository.save(pipelineOutput);

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");
    // mock GCS service
    when(mockGcsService.generateGetObjectSignedUrl(
            eq(pipelineRun.getWorkspaceGoogleProject()),
            eq(pipelineRun.getWorkspaceStorageContainerName()),
            anyString()))
        .thenReturn(fakeUrl);

    ApiPipelineRunOutputSignedUrls apiPipelineRunOutputs =
        pipelineInputsOutputsService.generatePipelineRunOutputSignedUrls(pipelineRun);

    assertEquals(fakeUrl.toString(), apiPipelineRunOutputs.get("testFileOutputKey"));
    // response should only include the file's signed url, not the other string output
    assertEquals(1, apiPipelineRunOutputs.size());
  }

  @Test
  void incrementDownloadCallCount() {
    UUID jobId = TestUtils.TEST_NEW_UUID;

    int updatedCallCount = pipelineInputsOutputsService.incrementDownloadCallCount(jobId);

    assertEquals(1, updatedCallCount);

    DownloadOutputCall retrievedCallEntity = downloadOutputCallRepository.findByJobId(jobId).get();
    assertEquals(1, retrievedCallEntity.getDownloadCallCount());
    assertEquals(
        retrievedCallEntity.getFirstCallTimestamp(), retrievedCallEntity.getLatestCallTimestamp());

    // now increment again to verify it updates existing record
    updatedCallCount = pipelineInputsOutputsService.incrementDownloadCallCount(jobId);
    assertEquals(2, updatedCallCount);

    retrievedCallEntity = downloadOutputCallRepository.findByJobId(jobId).get();
    assertEquals(2, retrievedCallEntity.getDownloadCallCount());
    assertTrue(
        retrievedCallEntity
            .getLatestCallTimestamp()
            .isAfter(retrievedCallEntity.getFirstCallTimestamp()));
  }

  @Test
  void stringToMapBadString() {
    assertThrows(
        InternalServerErrorException.class,
        () -> pipelineInputsOutputsService.stringToMap("i am not a map"));
  }

  @Test
  void getPipelineInputDefinitions() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    List<PipelineInputDefinition> allPipelineInputDefinitions =
        pipelinesService.getPipeline(imputationPipeline, null, false).getPipelineInputDefinitions();
    List<PipelineInputDefinition> userProvidedPipelineInputDefinitions =
        pipelineInputsOutputsService.extractUserProvidedInputDefinitions(
            allPipelineInputDefinitions);
    List<PipelineInputDefinition> serviceProvidedPipelineInputDefinitions =
        pipelineInputsOutputsService.extractServiceProvidedInputDefinitions(
            allPipelineInputDefinitions);

    assertEquals(
        allPipelineInputDefinitions.size(),
        userProvidedPipelineInputDefinitions.size()
            + serviceProvidedPipelineInputDefinitions.size());
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
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_STRING_INPUT_NAME,
                    "the-basename-value-for-my-output",
                    REQUIRED_VCF_INPUT_NAME,
                    "this/is/a/vcf/path.vcf.gz",
                    "extra_input_not_defined_in_pipeline",
                    "some value")),
            false,
            List.of(
                "Problem with pipelineInputs:",
                "Found extra input (extra_input_not_defined_in_pipeline)")),
        arguments(
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_STRING_INPUT_NAME,
                    "the-basename-value-for-my-output",
                    REQUIRED_VCF_INPUT_NAME,
                    "this/is/a/vcf/path.vcf.gz",
                    "extra_input_not_defined_in_pipeline",
                    "some value",
                    "another_extra_input",
                    "some_other_value")),
            false,
            List.of(
                "Problem with pipelineInputs:",
                "Found extra inputs",
                "extra_input_not_defined_in_pipeline",
                "another_extra_input")),
        arguments(
            new HashMap<String, Object>(Map.of()),
            false,
            List.of(
                "Problems with pipelineInputs:",
                "%s is required".formatted(REQUIRED_VCF_INPUT_NAME),
                "%s is required".formatted(REQUIRED_STRING_INPUT_NAME))),
        arguments(
            new HashMap<String, Object>(
                Map.of(REQUIRED_STRING_INPUT_NAME, List.of("this is an array, not a string"))),
            false,
            List.of(
                "Problems with pipelineInputs:",
                "%s is required".formatted(REQUIRED_VCF_INPUT_NAME),
                "%s must be a string".formatted(REQUIRED_STRING_INPUT_NAME))));
  }

  @ParameterizedTest
  @MethodSource("inputValidations")
  void validateInputs(
      Map<String, Object> inputs,
      Boolean shouldPassValidation,
      List<String> expectedErrorMessageStrings) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    List<PipelineInputDefinition> allInputDefinitions =
        pipelinesService.getPipeline(pipelinesEnum, null, false).getPipelineInputDefinitions();

    if (shouldPassValidation) {
      assertDoesNotThrow(
          () ->
              pipelineInputsOutputsService.validateUserProvidedInputs(allInputDefinitions, inputs));
    } else {
      ValidationException exception =
          assertThrows(
              ValidationException.class,
              () ->
                  pipelineInputsOutputsService.validateUserProvidedInputs(
                      allInputDefinitions, inputs));
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
            "input name",
            "description",
            PipelineVariableTypesEnum.INTEGER,
            null,
            isRequired,
            true,
            false,
            null,
            null,
            null);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    if (shouldPassValidation) {
      assertTrue(
          pipelineInputsOutputsService.validateRequiredInputs(inputDefinitions, inputs).isEmpty());
    } else {
      assertEquals(
          1, pipelineInputsOutputsService.validateRequiredInputs(inputDefinitions, inputs).size());
      assertEquals(
          "inputName is required",
          pipelineInputsOutputsService.validateRequiredInputs(inputDefinitions, inputs).get(0));
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
        arguments(PipelineVariableTypesEnum.STRING, null, "I am a string", false),
        arguments(PipelineVariableTypesEnum.STRING, null, "IAmAString", true),
        arguments(PipelineVariableTypesEnum.STRING, null, "  IAmAString  ", true),
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
            1L,
            "inputName",
            "input_name",
            "input name",
            "description",
            inputType,
            fileSuffix,
            true,
            true,
            false,
            null,
            null,
            null);
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>(List.of(inputDefinition));

    Map<String, Object> inputs = new HashMap<>();
    inputs.put("inputName", inputValue);

    // error message contents are tested in PipelineInputTypesEnumTest
    assertEquals(
        shouldPassValidation,
        pipelineInputsOutputsService.validateInputTypes(inputDefinitions, inputs).isEmpty());
  }

  private static Stream<Arguments> populateDefaultValuesForMissingOptionalUserInputsTestValues() {
    return Stream.of(
        // arguments: user-provided inputs, all pipeline input definitions, expected populated
        // inputs
        arguments( // one required user input, no optional/default values
            Map.of("inputName", "user provided value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, false, true, false, null)),
            Map.of("inputName", "user provided value")),
        arguments( // optional user input, not provided, uses default
            Map.of(),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, false, true, false, "default value")),
            Map.of("inputName", "default value")),
        arguments( // optional user input, provided, uses user value
            Map.of("inputName", "user provided value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, false, true, false, "default value")),
            Map.of("inputName", "user provided value")),
        arguments( // multiple input definitions
            Map.of(
                "inputName",
                "user provided value",
                "inputNameOptionalSpecified",
                "user specified value"),
            List.of(
                createTestPipelineInputDefWithName(
                    "inputName",
                    "input_name",
                    PipelineVariableTypesEnum.STRING,
                    false,
                    true,
                    false,
                    null,
                    null,
                    null),
                createTestPipelineInputDefWithName(
                    "inputNameOptionalNotSpecified",
                    "input_name_optional_not_specified",
                    PipelineVariableTypesEnum.STRING,
                    false,
                    true,
                    false,
                    "default optional value not specified by user",
                    null,
                    null),
                createTestPipelineInputDefWithName(
                    "inputNameOptionalSpecified",
                    "input_name_optional_specified",
                    PipelineVariableTypesEnum.STRING,
                    false,
                    true,
                    false,
                    "default optional value should be overridden by user value",
                    null,
                    null)),
            Map.of(
                "inputName",
                "user provided value",
                "inputNameOptionalNotSpecified",
                "default optional value not specified by user",
                "inputNameOptionalSpecified",
                "user specified value")));
  }

  @ParameterizedTest
  @MethodSource("populateDefaultValuesForMissingOptionalUserInputsTestValues")
  void populateDefaultValuesForMissingOptionalUserInputs(
      Map<String, Object> userProvidedInputs,
      List<PipelineInputDefinition> allPipelineInputDefinitions,
      Map<String, Object> expectedPopulatedInputs) {

    Map<String, Object> populatedInputs =
        pipelineInputsOutputsService.populateDefaultValuesForMissingOptionalUserInputs(
            allPipelineInputDefinitions, userProvidedInputs);

    assertEquals(expectedPopulatedInputs.size(), populatedInputs.size());
    for (String inputName :
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())) {
      assertEquals(expectedPopulatedInputs.get(inputName), populatedInputs.get(inputName));
    }
  }

  private static Stream<Arguments> addServiceProvidedInputsTestValues() {
    return Stream.of(
        // arguments: user-provided inputs, all pipeline input definitions, expected populated
        // inputs
        arguments( // one user input, no service provided inputs
            Map.of("inputName", "user provided value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, true, true, false, null)),
            Map.of("inputName", "user provided value")),
        arguments( // service provided input, use default value
            Map.of(),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING,
                    true,
                    false,
                    false,
                    "service provided default value")),
            Map.of("inputName", "service provided default value")),
        arguments( // service provided input with custom value
            Map.of(),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, true, false, true, null)),
            Collections.singletonMap("inputName", null)),
        arguments( // multiple input definitions
            Map.of("inputName", "user provided value"),
            List.of(
                createTestPipelineInputDefWithName(
                    "inputName",
                    "input_name",
                    PipelineVariableTypesEnum.STRING,
                    true,
                    true,
                    false,
                    null,
                    null,
                    null),
                createTestPipelineInputDefWithName(
                    "inputNameServiceProvided",
                    "input_name_service_provided",
                    PipelineVariableTypesEnum.STRING,
                    true,
                    false,
                    false,
                    "service provided default value",
                    null,
                    null)),
            Map.of(
                "inputName",
                "user provided value",
                "inputNameServiceProvided",
                "service provided default value")));
  }

  @ParameterizedTest
  @MethodSource("addServiceProvidedInputsTestValues")
  void addServiceProvidedInputs(
      Map<String, Object> userProvidedInputs,
      List<PipelineInputDefinition> allPipelineInputDefinitions,
      Map<String, Object> expectedPopulatedInputs) {

    Map<String, Object> populatedInputs =
        pipelineInputsOutputsService.addServiceProvidedInputs(
            allPipelineInputDefinitions, userProvidedInputs);

    assertEquals(expectedPopulatedInputs.size(), populatedInputs.size());
    for (String inputName :
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getName)
            .collect(Collectors.toSet())) {
      assertEquals(expectedPopulatedInputs.get(inputName), populatedInputs.get(inputName));
    }
  }

  private static Stream<Arguments> formatPipelineInputsTestValues() {
    return Stream.of(
        // arguments: all raw inputs, all pipeline input definitions, inputs with custom values,
        // keys to prepend with storage url, expected formatted inputs
        arguments( // no special modifications
            Map.of("inputName", "user provided value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, true, true, false, null)),
            Map.of(), // no inputs with custom values
            List.of(), // no keys to prepend with storage workspace url
            Map.of("input_name", "user provided value")),
        arguments( // overwrite service input value with custom value
            Collections.singletonMap("inputName", null),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, true, false, true, null)),
            Map.of("inputName", "custom value"),
            List.of(), // no keys to prepend with storage workspace url
            Map.of("input_name", "custom value")),
        arguments( // format user file with control workspace url and user input file path
            Map.of("inputName", "value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.FILE, true, true, false, null)),
            Map.of(), // no inputs with custom values
            List.of(), // no keys to prepend with storage workspace url
            Map.of(
                "input_name",
                "gs://control-workspace-bucket/user-input-files/%s/value"
                    .formatted(TestUtils.TEST_NEW_UUID))),
        arguments( // prepend key with storage workspace url
            Map.of("inputName", "/value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.STRING, true, false, false, null)),
            Map.of(), // no inputs with custom values
            List.of("inputName"), // prepend this key with storage workspace url
            Map.of("input_name", "gs://storage-workspace-bucket/value")),
        arguments( // test casting
            Map.of("inputName", "42"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.INTEGER, true, true, false, null)),
            Map.of(), // no inputs with custom values
            List.of(), // no keys to prepend with storage workspace url
            Map.of("input_name", 42)),
        arguments( // can handle multiple input definitions
            Map.of(
                "inputNameUserProvided",
                "user provided value",
                "inputNameServiceProvided",
                "service provided value"),
            List.of(
                createTestPipelineInputDefWithName(
                    "inputNameUserProvided",
                    "input_name_user_provided",
                    PipelineVariableTypesEnum.STRING,
                    true,
                    true,
                    false,
                    null,
                    null,
                    null),
                createTestPipelineInputDefWithName(
                    "inputNameServiceProvided",
                    "input_name_service_provided",
                    PipelineVariableTypesEnum.STRING,
                    true,
                    false,
                    false,
                    "service provided value",
                    null,
                    null)),
            Map.of(), // no inputs with custom values
            List.of(), // no keys to prepend with storage workspace url
            Map.of(
                "input_name_user_provided",
                "user provided value",
                "input_name_service_provided",
                "service provided value")));
  }

  @ParameterizedTest
  @MethodSource("formatPipelineInputsTestValues")
  void formatPipelineInputs(
      Map<String, Object> allRawInputs,
      List<PipelineInputDefinition> allPipelineInputDefinitions,
      Map<String, String> inputsWithCustomValues,
      List<String> keysToPrependWithStorageWorkspaceContainerUrl,
      Map<String, Object> expectedFormattedInputs) {

    UUID jobId = TestUtils.TEST_NEW_UUID;
    String controlWorkspaceContainerUrl = "gs://control-workspace-bucket";
    String storageWorkspaceContainerUrl = "gs://storage-workspace-bucket";

    Map<String, Object> formattedInputs =
        pipelineInputsOutputsService.formatPipelineInputs(
            allRawInputs,
            allPipelineInputDefinitions,
            jobId,
            controlWorkspaceContainerUrl,
            inputsWithCustomValues,
            keysToPrependWithStorageWorkspaceContainerUrl,
            storageWorkspaceContainerUrl);

    assertEquals(expectedFormattedInputs.size(), formattedInputs.size());

    for (String wdlVariableName :
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())) {
      assertEquals(
          expectedFormattedInputs.get(wdlVariableName), formattedInputs.get(wdlVariableName));
    }
  }

  @Test
  void gatherAndFormatPipelineInputs() {
    // test a couple bits of logic here, but the meat of the logic is tested separately above
    UUID jobId = TestUtils.TEST_NEW_UUID;
    String controlWorkspaceContainerUrl = "gs://control-workspace-bucket";
    String storageWorkspaceContainerUrl = "gs://storage-workspace-bucket";

    Map<String, Object> userInputs = Map.of("userInputName", "value");
    List<PipelineInputDefinition> allPipelineInputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "userInputName",
                "user_input_name",
                PipelineVariableTypesEnum.FILE,
                true,
                true,
                false,
                null,
                null,
                null),
            createTestPipelineInputDefWithName(
                "serviceInputName",
                "service_input_name",
                PipelineVariableTypesEnum.INTEGER,
                true,
                false,
                true,
                null,
                null,
                null));
    Map<String, String> inputsWithCustomValues = Map.of("serviceInputName", "123");
    List<String> keysToPrependWithStorageWorkspaceContainerUrl =
        List.of(); // no keys to prepend with storage workspace url

    Map<String, Object> expectedFormattedOutputs =
        Map.of(
            "user_input_name",
            "gs://control-workspace-bucket/user-input-files/%s/value".formatted(jobId),
            "service_input_name",
            123);

    Map<String, Object> formattedInputs =
        pipelineInputsOutputsService.gatherAndFormatPipelineInputs(
            jobId,
            allPipelineInputDefinitions,
            userInputs,
            controlWorkspaceContainerUrl,
            inputsWithCustomValues,
            keysToPrependWithStorageWorkspaceContainerUrl,
            storageWorkspaceContainerUrl);

    for (String wdlVariableName :
        allPipelineInputDefinitions.stream()
            .map(PipelineInputDefinition::getWdlVariableName)
            .collect(Collectors.toSet())) {
      assertEquals(
          expectedFormattedOutputs.get(wdlVariableName), formattedInputs.get(wdlVariableName));
    }
  }

  // test helper methods
  private static PipelineInputDefinition createTestPipelineInputDef(
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided,
      boolean isCustomValue,
      String defaultValue) {
    return createTestPipelineInputDefWithName(
        "inputName",
        "input_name",
        type,
        isRequired,
        isUserProvided,
        isCustomValue,
        defaultValue,
        null,
        null);
  }

  private static PipelineInputDefinition createTestPipelineInputDefWithName(
      String inputName,
      String inputWdlVariableName,
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided,
      boolean isCustomValue,
      String defaultValue,
      Double minValue,
      Double maxValue) {
    return new PipelineInputDefinition(
        3L,
        inputName,
        inputWdlVariableName,
        null,
        null,
        type,
        null,
        isRequired,
        isUserProvided,
        isCustomValue,
        defaultValue,
        minValue,
        maxValue);
  }
}
