package bio.terra.pipelines.service;

import static bio.terra.pipelines.testutils.TestUtils.createTestPipelineWithId;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyString;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.*;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.exception.ValidationException;
import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutput;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.gcs.GcsService;
import bio.terra.pipelines.dependencies.sam.SamService;
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
import org.mockito.ArgumentCaptor;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class PipelineInputsOutputsServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelineInputsOutputsService pipelineInputsOutputsService;

  @Autowired PipelinesService pipelinesService;
  @Autowired PipelineInputsRepository pipelineInputsRepository;
  @Autowired PipelineOutputsRepository pipelineOutputsRepository;

  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @MockitoBean private SamService mockSamService;
  @MockitoBean private GcsService mockGcsService;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final String fileInputKeyName1 = "testRequiredVcfInput";
  private final String fileInputKeyName2 = "testRequiredVcfInput2";
  private final List<PipelineInputDefinition> inputDefinitionsWithTwoFiles =
      List.of(
          new PipelineInputDefinition(
              3L,
              fileInputKeyName1,
              "test_required_vcf_input",
              null,
              null,
              PipelineVariableTypesEnum.FILE,
              ".vcf.gz",
              true,
              true,
              false,
              null,
              null,
              null),
          new PipelineInputDefinition(
              3L,
              fileInputKeyName2,
              "test_required_vcf_input_2",
              null,
              null,
              PipelineVariableTypesEnum.FILE,
              ".vcf.gz",
              true,
              true,
              false,
              null,
              null,
              null));

  @Test
  void validateFileSourcesAreConsistentCloudTrue() {
    String fileInputValue1 = "gs://some-bucket/some-path/file.vcf.gz";
    String fileInputValue2 = "gs://some-bucket/some-path/another-file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(
            Map.of(fileInputKeyName1, fileInputValue1, fileInputKeyName2, fileInputValue2));

    assertEquals(
        List.of(),
        pipelineInputsOutputsService.validateFileSourcesAreConsistent(
            inputDefinitionsWithTwoFiles, userPipelineInputs));
  }

  @Test
  void validateFileSourcesAreConsistentLocalTrue() {
    String fileInputValue1 = "some-path/file.vcf.gz";
    String fileInputValue2 = "some-path/another-file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(
            Map.of(fileInputKeyName1, fileInputValue1, fileInputKeyName2, fileInputValue2));

    assertEquals(
        List.of(),
        pipelineInputsOutputsService.validateFileSourcesAreConsistent(
            inputDefinitionsWithTwoFiles, userPipelineInputs));
  }

  @Test
  void validateFileSourcesAreConsistentMixError() {
    String fileInputValue1 = "gs://some-bucket/some-path/file.vcf.gz";
    String fileInputValue2 = "/some-path/another-file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(
            Map.of(fileInputKeyName1, fileInputValue1, fileInputKeyName2, fileInputValue2));

    assertEquals(
        List.of("File inputs must be all local or all GCS cloud based"),
        pipelineInputsOutputsService.validateFileSourcesAreConsistent(
            inputDefinitionsWithTwoFiles, userPipelineInputs));
  }

  @Test
  void validateFileSourcesAreConsistentNonGcsCloudErrors() {
    String fileInputValue1 = "s3://some-bucket/some-path/file.vcf.gz";
    String fileInputValue2 = "azure://some-path/another-file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(
            Map.of(fileInputKeyName1, fileInputValue1, fileInputKeyName2, fileInputValue2));

    assertEquals(
        List.of(
            "Found an unsupported file location type for input %s. Only GCS cloud-based files or local files are supported"
                .formatted(fileInputKeyName1),
            "Found an unsupported file location type for input %s. Only GCS cloud-based files or local files are supported"
                .formatted(fileInputKeyName2)),
        pipelineInputsOutputsService.validateFileSourcesAreConsistent(
            inputDefinitionsWithTwoFiles, userPipelineInputs));
  }

  private static final String FILE_WITH_USER_ACCESS_ONLY =
      "gs://bucket/file-with-user-access-only.vcf.gz";
  private static final String FILE_WITH_SERVICE_ACCESS_ONLY =
      "gs://bucket/file-with-service-access-only.vcf.gz";
  private static final String FILE_WITH_SERVICE_AND_USER_ACCESS =
      "gs://bucket/file-with-both-access.vcf.gz";
  private static final String FILE_WITH_NO_ACCESS = "gs://bucket/file-with-no-access.vcf.gz";
  private static final String USER_PROXY_GROUP = "PROXY_pizza@cake.com";

  private static final String USER_ACCESS_ERROR_MESSAGE_FORMAT =
      "User does not have necessary permissions to access file input for %s (%s), or the file does not exist. Please ensure the user's proxy group %s has read access to the bucket containing all input files, or that the files exist if the permissions are correct.";
  private static final String SERVICE_ACCESS_ERROR_MESSAGE_FORMAT =
      "Service does not have necessary permissions to access file input for %s (%s), or the file does not exist. Please ensure that test-service-account-group@test.com has read access to the bucket containing all input files, or that the files exist if the permissions are correct.";

  private static Stream<Arguments> validateAccessInputs() {
    return Stream.of(
        // arguments: userProvidedInputs, shouldPassValidation, expectedErrorMessageString
        arguments( // user and service have access
            new HashMap<String, Object>(
                Map.of(
                    "testRequiredVcfInput",
                    FILE_WITH_SERVICE_AND_USER_ACCESS,
                    "testOptionalVcfInput",
                    FILE_WITH_SERVICE_AND_USER_ACCESS)),
            true,
            null),
        arguments( // user and service have access, missing an input is ok here
            new HashMap<String, Object>(
                Map.of("testRequiredVcfInput", FILE_WITH_SERVICE_AND_USER_ACCESS)),
            true,
            null),
        arguments( // user does not have access
            new HashMap<String, Object>(
                Map.of("testRequiredVcfInput", FILE_WITH_SERVICE_ACCESS_ONLY)),
            false,
            USER_ACCESS_ERROR_MESSAGE_FORMAT.formatted(
                "testRequiredVcfInput", FILE_WITH_SERVICE_ACCESS_ONLY, USER_PROXY_GROUP)),
        arguments( // service does not have access
            new HashMap<String, Object>(Map.of("testRequiredVcfInput", FILE_WITH_USER_ACCESS_ONLY)),
            false,
            SERVICE_ACCESS_ERROR_MESSAGE_FORMAT.formatted(
                "testRequiredVcfInput", FILE_WITH_USER_ACCESS_ONLY)),
        arguments( // neither user nor service have access; user check happens first
            new HashMap<String, Object>(Map.of("testRequiredVcfInput", FILE_WITH_NO_ACCESS)),
            false,
            USER_ACCESS_ERROR_MESSAGE_FORMAT.formatted(
                "testRequiredVcfInput", FILE_WITH_NO_ACCESS, USER_PROXY_GROUP)),
        arguments( // one input file has ok access but the other doesn't - should still fail
            new HashMap<String, Object>(
                Map.of(
                    "testRequiredVcfInput",
                    FILE_WITH_SERVICE_AND_USER_ACCESS,
                    "testOptionalVcfInput",
                    FILE_WITH_NO_ACCESS)),
            false,
            USER_ACCESS_ERROR_MESSAGE_FORMAT.formatted(
                "testOptionalVcfInput", FILE_WITH_NO_ACCESS, USER_PROXY_GROUP)),
        arguments( // not a GCS file logs an error but doesn't throw
            new HashMap<String, Object>(Map.of("testRequiredVcfInput", "not-a-gcs-file.vcf.gz")),
            true,
            null));
  }

  @ParameterizedTest
  @MethodSource("validateAccessInputs")
  void validateUserAndServiceReadAccessToCloudInputs(
      Map<String, Object> userProvidedInputs,
      boolean shouldPassValidation,
      String expectedErrorMessageString) {
    Pipeline pipeline = createTestPipelineWithId();
    // create inputs with one required file input, one optional, and one service-provided file input
    pipeline.setPipelineInputDefinitions(
        List.of(
            createTestPipelineInputDefWithName(
                "testRequiredVcfInput",
                "test_required_vcf_input",
                PipelineVariableTypesEnum.FILE,
                true,
                true),
            createTestPipelineInputDefWithName(
                "testOptionalVcfInput",
                "test_optional_vcf_input",
                PipelineVariableTypesEnum.FILE,
                false,
                true),
            createTestPipelineInputDefWithName(
                "testServiceProvidedVcfInput",
                "test_service_provided_vcf_input",
                PipelineVariableTypesEnum.FILE,
                true,
                false)));

    SamUser testUser = TestUtils.TEST_SAM_USER_1;
    BearerToken userPetToken = testUser.getBearerToken();
    String userPetTokenString = userPetToken.getToken();
    // mock the GCS service to return the expected access results
    when(mockGcsService.userHasFileReadAccess(
            new GcsFile(FILE_WITH_USER_ACCESS_ONLY), userPetTokenString))
        .thenReturn(true);
    when(mockGcsService.userHasFileReadAccess(
            new GcsFile(FILE_WITH_SERVICE_AND_USER_ACCESS), userPetTokenString))
        .thenReturn(true);
    when(mockGcsService.userHasFileReadAccess(
            new GcsFile(FILE_WITH_SERVICE_ACCESS_ONLY), userPetTokenString))
        .thenReturn(false);
    when(mockGcsService.userHasFileReadAccess(new GcsFile(FILE_WITH_NO_ACCESS), userPetTokenString))
        .thenReturn(false);

    when(mockGcsService.serviceHasFileReadAccess(new GcsFile(FILE_WITH_SERVICE_ACCESS_ONLY)))
        .thenReturn(true);
    when(mockGcsService.serviceHasFileReadAccess(new GcsFile(FILE_WITH_SERVICE_AND_USER_ACCESS)))
        .thenReturn(true);
    when(mockGcsService.serviceHasFileReadAccess(new GcsFile(FILE_WITH_USER_ACCESS_ONLY)))
        .thenReturn(false);
    when(mockGcsService.serviceHasFileReadAccess(new GcsFile(FILE_WITH_NO_ACCESS)))
        .thenReturn(false);

    when(mockSamService.getProxyGroupForUser(testUser)).thenReturn(USER_PROXY_GROUP);
    when(mockSamService.getUserPetServiceAccountTokenReadOnly(testUser)).thenReturn(userPetToken);

    if (shouldPassValidation) {
      assertDoesNotThrow(
          () ->
              pipelineInputsOutputsService.validateUserAndServiceReadAccessToCloudInputs(
                  pipeline, userProvidedInputs, testUser));
    } else {
      ValidationException exception =
          assertThrows(
              ValidationException.class,
              () ->
                  pipelineInputsOutputsService.validateUserAndServiceReadAccessToCloudInputs(
                      pipeline, userProvidedInputs, testUser));
      assertTrue(exception.getMessage().contains(expectedErrorMessageString));
    }
  }

  @Test
  void prepareLocalFileInputs() throws MalformedURLException {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generatePutObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceStorageContainerName()), anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineInputsOutputsService.prepareLocalFileInputs(
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
  void prepareLocalFileInputsResumable() throws MalformedURLException {
    Pipeline testPipelineWithId = createTestPipelineWithId();
    String fileInputKeyName = "testRequiredVcfInput";
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(fileInputKeyName, fileInputValue));

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generateResumablePostObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceStorageContainerName()), anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineInputsOutputsService.prepareLocalFileInputs(
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

    URL fakeSignedUrl = new URL("https://storage.googleapis.com/signed-url-stuff");
    // mock GCS service
    when(mockGcsService.generateGetObjectSignedUrl(any(GcsFile.class))).thenReturn(fakeSignedUrl);

    ApiPipelineRunOutputSignedUrls apiPipelineRunOutputs =
        pipelineInputsOutputsService.generatePipelineRunOutputSignedUrls(pipelineRun);

    assertEquals(fakeSignedUrl.toString(), apiPipelineRunOutputs.get("testFileOutputKey"));
    // response should only include the file's signed url, not the other string output
    assertEquals(1, apiPipelineRunOutputs.size());
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
  static final String OPTIONAL_VCF_INPUT_NAME = "optionalVcfInput";

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
                "%s must be a string".formatted(REQUIRED_STRING_INPUT_NAME))),
        arguments(
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_VCF_INPUT_NAME,
                    "this/is/a/vcf/path.vcf.gz", // local path
                    OPTIONAL_VCF_INPUT_NAME,
                    "gs://some-bucket/some-path/file.vcf.gz" // cloud path
                    )),
            false,
            List.of(
                "Problems with pipelineInputs:",
                "File inputs must be all local or all GCS cloud based")),
        arguments(
            new HashMap<String, Object>(
                Map.of(
                    REQUIRED_VCF_INPUT_NAME,
                    "this/is/a/vcf/path.vcf.gz", // local path
                    OPTIONAL_VCF_INPUT_NAME,
                    "s3://some-bucket/some-path/file.vcf.gz" // cloud path
                    )),
            false,
            List.of(
                "Problems with pipelineInputs:",
                "Found an unsupported file location type for input %s. Only GCS cloud-based files or local files are supported"
                    .formatted(OPTIONAL_VCF_INPUT_NAME))));
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

    // add optional file input for testing mixed local/cloud validation
    allInputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            OPTIONAL_VCF_INPUT_NAME,
            "optional_input_vcf",
            "optional input vcf",
            "description",
            PipelineVariableTypesEnum.FILE,
            ".vcf.gz",
            false,
            true,
            false,
            null,
            null,
            null));

    if (shouldPassValidation) {
      assertDoesNotThrow(
          () ->
              pipelineInputsOutputsService.validateUserProvidedInputsWithCloud(
                  allInputDefinitions, inputs));
    } else {
      ValidationException exception =
          assertThrows(
              ValidationException.class,
              () ->
                  pipelineInputsOutputsService.validateUserProvidedInputsWithCloud(
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
        arguments( // don't format cloud-based user file with control workspace url
            Map.of("inputName", "gs://bucket/value"),
            List.of(
                createTestPipelineInputDef(
                    PipelineVariableTypesEnum.FILE, true, true, false, null)),
            Map.of(), // no inputs with custom values
            List.of(), // no keys to prepend with storage workspace url
            Map.of("input_name", "gs://bucket/value")),
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
      boolean isUserProvided) {
    return createTestPipelineInputDefWithName(
        inputName, inputWdlVariableName, type, isRequired, isUserProvided, false, null, null, null);
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
    String fileSuffix =
        type == PipelineVariableTypesEnum.FILE || type == PipelineVariableTypesEnum.FILE_ARRAY
            ? ".vcf.gz"
            : null;
    return new PipelineInputDefinition(
        3L,
        inputName,
        inputWdlVariableName,
        null,
        null,
        type,
        fileSuffix,
        isRequired,
        isUserProvided,
        isCustomValue,
        defaultValue,
        minValue,
        maxValue);
  }

  // deliverOutputFilesToGcs tests

  @Test
  void deliverOutputFilesToGcsSuccess() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);
    GcsFile destinationGcsPath = new GcsFile("gs://destination-bucket/path");

    // Create test outputs
    Map<String, Object> outputsMap = new HashMap<>();
    outputsMap.put("outputFile", "gs://source-bucket/path/to/file.vcf.gz");

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    ArgumentCaptor<GcsFile> sourceFileCaptor = ArgumentCaptor.forClass(GcsFile.class);
    ArgumentCaptor<GcsFile> destinationFileCaptor = ArgumentCaptor.forClass(GcsFile.class);

    doNothing().when(mockGcsService).copyObject(any(GcsFile.class), any(GcsFile.class));

    pipelineInputsOutputsService.deliverOutputFilesToGcs(testPipelineRun, destinationGcsPath);

    // Verify copyObject was called and capture the arguments
    verify(mockGcsService).copyObject(sourceFileCaptor.capture(), destinationFileCaptor.capture());

    // Verify the source file
    GcsFile capturedSourceFile = sourceFileCaptor.getValue();
    assertEquals("gs://source-bucket/path/to/file.vcf.gz", capturedSourceFile.getFullPath());

    // Verify the destination path includes the jobId folder
    GcsFile capturedDestinationFile = destinationFileCaptor.getValue();
    String capturedDestinationPath = capturedDestinationFile.getFullPath();
    assertTrue(capturedDestinationPath.contains(testPipelineRun.getJobId().toString() + "/"));
    assertTrue(capturedDestinationPath.endsWith("file.vcf.gz"));
  }

  @Test
  void deliverOutputFilesToGcsWithJobIdFolder() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);
    GcsFile destinationGcsPath = new GcsFile("gs://destination-bucket/path");

    // Create test outputs
    Map<String, Object> outputsMap = new HashMap<>();
    outputsMap.put("outputFile", "gs://source-bucket/path/to/file.vcf.gz");

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    ArgumentCaptor<GcsFile> sourceFileCaptor = ArgumentCaptor.forClass(GcsFile.class);
    ArgumentCaptor<GcsFile> destinationFileCaptor = ArgumentCaptor.forClass(GcsFile.class);

    Mockito.doNothing().when(mockGcsService).copyObject(any(GcsFile.class), any(GcsFile.class));

    pipelineInputsOutputsService.deliverOutputFilesToGcs(testPipelineRun, destinationGcsPath);

    // Verify copyObject was called and capture the arguments
    Mockito.verify(mockGcsService)
        .copyObject(sourceFileCaptor.capture(), destinationFileCaptor.capture());

    // Verify the destination path includes the jobId folder
    GcsFile capturedDestinationFile = destinationFileCaptor.getValue();
    String capturedDestinationPath = capturedDestinationFile.getFullPath();
    assertTrue(capturedDestinationPath.contains(testPipelineRun.getJobId().toString() + "/"));
    assertTrue(capturedDestinationPath.endsWith("file.vcf.gz"));
  }

  @Test
  void deliverOutputFilesToGcsCopyFailureThrowsException() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);
    GcsFile destinationGcsPath = new GcsFile("gs://destination-bucket/path");

    // Create test output
    Map<String, Object> outputsMap = new HashMap<>();
    outputsMap.put("outputFile", "gs://source-bucket/path/to/file.vcf.gz");

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    // Mock GCS service to throw an exception
    doThrow(new RuntimeException("oh no something broke!"))
        .when(mockGcsService)
        .copyObject(any(GcsFile.class), any(GcsFile.class));

    // Should throw InternalServerErrorException
    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineInputsOutputsService.deliverOutputFilesToGcs(
                testPipelineRun, destinationGcsPath));
  }

  @Test
  void deliverOutputFilesToGcsWithEmptyOutputsMap() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);
    GcsFile destinationGcsPath = new GcsFile("gs://destination-bucket/path");

    // Create empty outputs
    Map<String, Object> outputsMap = new HashMap<>();

    // Save outputs to repository
    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    // Should not throw an exception even with empty outputs
    assertDoesNotThrow(
        () ->
            pipelineInputsOutputsService.deliverOutputFilesToGcs(
                testPipelineRun, destinationGcsPath));

    // Verify that copyObject was never called
    verify(mockGcsService, never()).copyObject(any(GcsFile.class), any(GcsFile.class));
  }

  // Helper method to create a test pipeline run
  private PipelineRun createTestPipelineRun(Pipeline pipeline, UUID jobId) {
    PipelineRun pipelineRun = new PipelineRun();
    pipelineRun.setJobId(jobId);
    pipelineRun.setUserId("test-user");
    pipelineRun.setPipelineId(pipeline.getId());
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRun.setWorkspaceBillingProject(pipeline.getWorkspaceBillingProject());
    pipelineRun.setWorkspaceName(pipeline.getWorkspaceName());
    pipelineRun.setWorkspaceStorageContainerName(pipeline.getWorkspaceStorageContainerName());
    pipelineRun.setWorkspaceGoogleProject(pipeline.getWorkspaceGoogleProject());
    return pipelineRunsRepository.save(pipelineRun);
  }
}
