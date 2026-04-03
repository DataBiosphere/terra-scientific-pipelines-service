package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.FileUtils.constructDestinationBlobNameForUserInputFile;
import static bio.terra.pipelines.testutils.TestUtils.*;
import static bio.terra.pipelines.testutils.TestUtils.updateTestPipeline1WithTestValues;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.exception.ValidationException;
import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.*;
import bio.terra.pipelines.db.repositories.PipelineInputDefinitionsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputDefinitionsRepository;
import bio.terra.pipelines.db.repositories.PipelineOutputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
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
import java.util.Set;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.jupiter.api.BeforeEach;
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
  @Autowired PipelinesRepository pipelinesRepository;
  @Autowired PipelineInputDefinitionsRepository pipelineInputDefinitionsRepository;
  @Autowired PipelineOutputDefinitionsRepository pipelineOutputDefinitionsRepository;

  @MockitoBean private SamService mockSamService;
  @MockitoBean private GcsService mockGcsService;

  private static final UUID TEST_JOB_ID = TestUtils.TEST_NEW_UUID;
  private static final String FILE_INPUT_KEY_NAME = "testRequiredVcfInput";
  private static final String MANIFEST_INPUT_KEY_NAME = "testRequiredManifestInput";
  private static final List<PipelineInputDefinition> INPUT_DEFINITIONS_WITH_FILE_AND_MANIFEST =
      List.of(
          createTestPipelineInputDefWithName(
              FILE_INPUT_KEY_NAME,
              "test_required_vcf_input",
              PipelineVariableTypesEnum.FILE,
              true,
              true),
          createTestPipelineInputDefWithName(
              MANIFEST_INPUT_KEY_NAME,
              "test_required_manifest_input",
              PipelineVariableTypesEnum.MANIFEST,
              true,
              true));

  @BeforeEach
  void addTestPipelineToDb() {
    // add a test pipeline with one required file input and one required manifest input to the db,
    // so
    // we can use it for validation tests
    Pipeline testPipeline = addNewTestPipelineWithTestValues();
    Pipeline savedPipeline = pipelinesRepository.save(testPipeline);
    Long pipelineId = savedPipeline.getId();

    pipelineInputDefinitionsRepository.saveAll(
        List.of(
            createTestPipelineInputDefWithNameAndPipelineId(
                pipelineId,
                FILE_INPUT_KEY_NAME,
                "test_required_vcf_input",
                PipelineVariableTypesEnum.FILE,
                true,
                true),
            createTestPipelineInputDefWithNameAndPipelineId(
                pipelineId,
                MANIFEST_INPUT_KEY_NAME,
                "test_required_manifest_input",
                PipelineVariableTypesEnum.MANIFEST,
                true,
                true)));

    pipelineOutputDefinitionsRepository.saveAll(
        List.of(
            new PipelineOutputDefinition(
                pipelineId,
                "testFileOutputKey",
                "test_file_output_key",
                "Test File Output Display Name",
                "test output file description",
                PipelineVariableTypesEnum.FILE,
                true),
            new PipelineOutputDefinition(
                pipelineId,
                "testStringOutputKey",
                "test_string_output_key",
                "Test String Output Display Name",
                "test output string description",
                PipelineVariableTypesEnum.STRING,
                true)));
  }

  private static Stream<Arguments> validateFileSourcesAreConsistentTestInputs() {
    return Stream.of(
        // arguments: userProvidedInputs, expectedErrorMessageStrings (passes validation = empty
        // list)
        arguments( // two GCS cloud inputs
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/manifest.tsv")),
            List.of()),
        arguments( // two local inputs
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "some-path/manifest.tsv")),
            List.of()),
        arguments( // local and cloud
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "some-path/manifest.tsv")),
            List.of("File inputs must be all local or all GCS cloud based")),
        arguments( // mixed wrong cloud
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "s3://some-bucket/some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "azure://some-path/manifest.tsv")),
            List.of(
                "Found an unsupported file location type for input %s. Only GCS cloud-based files or local files are supported"
                    .formatted(FILE_INPUT_KEY_NAME),
                "Found an unsupported file location type for input %s. Only GCS cloud-based files or local files are supported"
                    .formatted(MANIFEST_INPUT_KEY_NAME))));
  }

  @ParameterizedTest
  @MethodSource("validateFileSourcesAreConsistentTestInputs")
  void validateFileSourcesAreConsistent(
      Map<String, Object> userPipelineInputs, List<String> expectedErrorMessageStrings) {
    assertEquals(
        expectedErrorMessageStrings,
        pipelineInputsOutputsService.validateFileSourcesAreConsistent(
            INPUT_DEFINITIONS_WITH_FILE_AND_MANIFEST, userPipelineInputs));
  }

  private static Stream<Arguments> userFileInputsAreCloudTestInputs() {
    return Stream.of(
        // arguments: userProvidedInputs, areCloud
        arguments( // two GCS cloud inputs
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/manifest.tsv")),
            true),
        arguments( // one local, one cloud
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "/some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/manifest.tsv")),
            false),
        arguments( // both local
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "manifest.tsv")),
            false),
        arguments( // GCS and non-GCS cloud
            new HashMap<>(
                Map.of(
                    FILE_INPUT_KEY_NAME,
                    "gs://some-bucket/some-path/file.vcf.gz",
                    MANIFEST_INPUT_KEY_NAME,
                    "s3://some-bucket/some-path/manifest.tsv")),
            false));
  }

  @ParameterizedTest
  @MethodSource("userFileInputsAreCloudTestInputs")
  void userProvidedInputsAreGcsCloud(Map<String, Object> userPipelineInputs, boolean areCloud) {
    Pipeline testPipeline = updateTestPipeline1WithTestValues();
    testPipeline.setPipelineInputDefinitions(INPUT_DEFINITIONS_WITH_FILE_AND_MANIFEST);

    if (areCloud) {
      assertTrue(
          pipelineInputsOutputsService.userProvidedInputsAreGcsCloud(
              testPipeline, userPipelineInputs));
    } else {
      assertFalse(
          pipelineInputsOutputsService.userProvidedInputsAreGcsCloud(
              testPipeline, userPipelineInputs));
    }
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
            null),
        arguments( // file and manifest with access
            new HashMap<String, Object>(
                Map.of(
                    "testRequiredVcfInput",
                    FILE_WITH_SERVICE_AND_USER_ACCESS,
                    "testRequiredManifestInput",
                    FILE_WITH_SERVICE_AND_USER_ACCESS)),
            true,
            null),
        arguments( // manifest file with user access but no service access
            new HashMap<String, Object>(
                Map.of(
                    "testRequiredVcfInput",
                    FILE_WITH_SERVICE_AND_USER_ACCESS,
                    "testRequiredManifestInput",
                    FILE_WITH_USER_ACCESS_ONLY)),
            false,
            SERVICE_ACCESS_ERROR_MESSAGE_FORMAT.formatted(
                "testRequiredManifestInput", FILE_WITH_USER_ACCESS_ONLY)));
  }

  @ParameterizedTest
  @MethodSource("validateAccessInputs")
  void validateUserAndServiceReadAccessToCloudInputs(
      Map<String, Object> userProvidedInputs,
      boolean shouldPassValidation,
      String expectedErrorMessageString) {
    Pipeline pipeline = updateTestPipeline1WithTestValues();
    // create inputs with one required file input, one required manifest input, one optional file
    // input, and one service-provided file input
    pipeline.setPipelineInputDefinitions(
        List.of(
            createTestPipelineInputDefWithName(
                "testRequiredVcfInput",
                "test_required_vcf_input",
                PipelineVariableTypesEnum.FILE,
                true,
                true),
            createTestPipelineInputDefWithName(
                "testRequiredManifestInput",
                "test_required_manifest_input",
                PipelineVariableTypesEnum.MANIFEST,
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
    Pipeline testPipelineWithId = updateTestPipeline1WithTestValues();
    testPipelineWithId.setPipelineInputDefinitions(INPUT_DEFINITIONS_WITH_FILE_AND_MANIFEST);
    String fileInputValue = "fake/file.vcf.gz";
    String manifestInputValue = "fake/manifest.tsv";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(
            Map.of(
                FILE_INPUT_KEY_NAME, fileInputValue, MANIFEST_INPUT_KEY_NAME, manifestInputValue));

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generatePutObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceStorageContainerName()), anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineInputsOutputsService.prepareLocalFileInputs(
            testPipelineWithId, TEST_JOB_ID, userPipelineInputs, false);

    assertEquals(userPipelineInputs.size(), formattedPipelineFileInputs.size());
    assertEquals(
        fakeUrl.toString(), formattedPipelineFileInputs.get(FILE_INPUT_KEY_NAME).get("signedUrl"));
    assertEquals(
        "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s' | cat"
            .formatted(fileInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(FILE_INPUT_KEY_NAME).get("curlCommand"));

    assertEquals(
        fakeUrl.toString(),
        formattedPipelineFileInputs.get(MANIFEST_INPUT_KEY_NAME).get("signedUrl"));
    assertEquals(
        "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s '%s' | cat"
            .formatted(manifestInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(MANIFEST_INPUT_KEY_NAME).get("curlCommand"));
  }

  @Test
  void prepareLocalFileInputsResumable() throws MalformedURLException {
    Pipeline testPipelineWithId = updateTestPipeline1WithTestValues();
    String fileInputValue = "fake/file.vcf.gz";
    Map<String, Object> userPipelineInputs =
        new HashMap<>(Map.of(FILE_INPUT_KEY_NAME, fileInputValue));

    URL fakeUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

    when(mockGcsService.generateResumablePostObjectSignedUrl(
            eq(testPipelineWithId.getWorkspaceStorageContainerName()), anyString()))
        .thenReturn(fakeUrl);

    Map<String, Map<String, String>> formattedPipelineFileInputs =
        pipelineInputsOutputsService.prepareLocalFileInputs(
            testPipelineWithId, TEST_JOB_ID, userPipelineInputs, true);

    assertEquals(userPipelineInputs.size(), formattedPipelineFileInputs.size());
    assertEquals(
        fakeUrl.toString(), formattedPipelineFileInputs.get(FILE_INPUT_KEY_NAME).get("signedUrl"));
    assertEquals(
        "curl --progress-bar -X PUT -H 'Content-Type: application/octet-stream' --upload-file %s $(curl -s -i -X POST -H 'x-goog-resumable: start' '%s' | grep -i '^Location:' | cut -d' ' -f2- | tr -d '\r') | cat"
            .formatted(fileInputValue, fakeUrl.toString()),
        formattedPipelineFileInputs.get(FILE_INPUT_KEY_NAME).get("curlCommand"));
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
  void getPipelineRunOutputsV2() {
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(TEST_JOB_ID);
    pipelineRun.setPipelineId(
        pipelinesService
            .getPipeline(TEST_PIPELINE_1_IMPUTATION_ENUM, TEST_PIPELINE_VERSION_1, false)
            .getId());
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRunsRepository.save(pipelineRun);

    pipelineOutputsRepository.saveAll(getPipelineOutputsForPipelineRun(pipelineRun, true));

    Map<String, Object> retrievedOutputs =
        pipelineInputsOutputsService.getPipelineRunOutputsV2(pipelineRun);

    assertEquals(TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.size(), retrievedOutputs.size());
    for (String outputKey : TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.keySet()) {
      assertEquals(
          TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED.get(outputKey),
          retrievedOutputs.get(outputKey));
    }
  }

  @Test
  void getPipelineRunOutputsV3() {
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(TEST_JOB_ID);
    pipelineRun.setPipelineId(
        pipelinesService
            .getPipeline(TEST_PIPELINE_1_IMPUTATION_ENUM, TEST_PIPELINE_VERSION_1, false)
            .getId());
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRunsRepository.save(pipelineRun);

    pipelineOutputsRepository.saveAll(getPipelineOutputsForPipelineRun(pipelineRun, true));

    Map<String, Object> retrievedOutputs =
        pipelineInputsOutputsService.getPipelineRunOutputsV3(pipelineRun);

    assertEquals(TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.size(), retrievedOutputs.size());

    // extract and assert file output
    Object fileOutputWithMetadata = retrievedOutputs.get("testFileOutputKey");
    assertTrue(fileOutputWithMetadata instanceof Map);
    Map<String, Object> fileOutputMap = (Map<String, Object>) fileOutputWithMetadata;

    assertEquals(
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED.get("testFileOutputKey"),
        fileOutputMap.get("value"));
    assertTrue(fileOutputMap.containsKey("metadata"));
    Map<String, Object> fileMetadata = (Map<String, Object>) fileOutputMap.get("metadata");
    assertEquals(256L, fileMetadata.get("sizeInBytes"));

    // extract and assert string output (no metadata)
    Map<String, String> expectedStringOutput =
        Map.of(
            "value",
            TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED.get("testStringOutputKey"));
    assertEquals(expectedStringOutput, retrievedOutputs.get("testStringOutputKey"));

    // assert string output doesn't contain metadata key
    Map<String, Object> stringOutputMap =
        (Map<String, Object>) retrievedOutputs.get("testStringOutputKey");
    assertFalse(stringOutputMap.containsKey("metadata"));
  }

  @Test
  void getPipelineRunOutputsV3WithoutFileSize() {
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(TEST_JOB_ID);
    pipelineRun.setPipelineId(
        pipelinesService
            .getPipeline(TEST_PIPELINE_1_IMPUTATION_ENUM, TEST_PIPELINE_VERSION_1, false)
            .getId());
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRunsRepository.save(pipelineRun);

    pipelineOutputsRepository.saveAll(getPipelineOutputsForPipelineRun(pipelineRun, false));

    Map<String, Object> retrievedOutputs =
        pipelineInputsOutputsService.getPipelineRunOutputsV3(pipelineRun);

    assertEquals(TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.size(), retrievedOutputs.size());

    // extract and assert file output
    Object fileOutputWithMetadata = retrievedOutputs.get("testFileOutputKey");
    assertTrue(fileOutputWithMetadata instanceof Map);
    Map<String, Object> fileOutputMap = (Map<String, Object>) fileOutputWithMetadata;

    assertEquals(
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED.get("testFileOutputKey"),
        fileOutputMap.get("value"));

    // if file size is not available, metadata should not be included in response
    assertFalse(fileOutputMap.containsKey("metadata"));

    // extract and assert string output (no metadata)
    Map<String, String> expectedStringOutput =
        Map.of(
            "value",
            TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED.get("testStringOutputKey"));
    assertEquals(expectedStringOutput, retrievedOutputs.get("testStringOutputKey"));

    // assert string output doesn't contain metadata key
    Map<String, Object> stringOutputMap =
        (Map<String, Object>) retrievedOutputs.get("testStringOutputKey");
    assertFalse(stringOutputMap.containsKey("metadata"));
  }

  @Test
  void generatePipelineRunOutputSignedUrls() throws MalformedURLException {
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(TEST_JOB_ID);
    pipelineRun.setPipelineId(
        pipelinesService
            .getPipeline(TEST_PIPELINE_1_IMPUTATION_ENUM, TEST_PIPELINE_VERSION_1, false)
            .getId());
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRunsRepository.save(pipelineRun);

    pipelineOutputsRepository.saveAll(getPipelineOutputsForPipelineRun(pipelineRun, true));

    URL fakeSignedUrl = new URL("https://storage.googleapis.com/signed-url-stuff");

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

  @Test
  void extractUniqueBucketsFromManifests() {
    // test multiple inputs, multiple manifests, don't act on FILE input
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true),
            createTestPipelineInputDefWithName(
                "manifest2", "manifest_2", PipelineVariableTypesEnum.MANIFEST, false, true),
            createTestPipelineInputDefWithName(
                "file1", "file_1", PipelineVariableTypesEnum.FILE, false, true));

    String manifestFile1 = "gs://bucket1/path/to/manifest1.tsv";
    String manifestFile2 = "gs://bucket2/path/to/manifest2.tsv";

    String file1 = "gs://bucket4/path/to/file1.vcf.gz";
    String file2 = "gs://bucket5/path/to/file2.vcf.gz";
    String file3 = "gs://bucket6/path/to/file3.vcf.gz";
    String file4 = "gs://bucket6/path/to/file4.vcf.gz";
    String file5 = "gs://bucket7/path/to/file5.vcf.gz";
    String file6 = "gs://bucket8/path/to/file6.vcf.gz";

    String manifestFile1Contents =
        "sample1\t%s\t%s%nsample2\t%s\t%s%n".formatted(file1, file2, file3, file4);
    String manifestFile2Contents = "sample3\t%s%nsample4\t%s".formatted(file5, file6);

    // expected buckets are the buckets from the files, not from the manifests
    Set<String> expectedBucketSet = Set.of("bucket4", "bucket5", "bucket6", "bucket7", "bucket8");

    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile1Contents));
    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile2)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile2Contents));

    Map<String, Object> userInputs =
        Map.of(
            "manifest1", manifestFile1,
            "manifest2", manifestFile2,
            "file1", "gs://bucket3/path/to/file.vcf.gz");

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    Set<String> uniqueBucketsResult =
        pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
            inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun);

    assertEquals(expectedBucketSet.size(), uniqueBucketsResult.size());
    assertEquals(expectedBucketSet, uniqueBucketsResult);
  }

  @Test
  void extractUniqueBucketsFromManifestsNoManifestsOk() {
    // no manifest input
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "file1", "file_1", PipelineVariableTypesEnum.FILE, false, true));

    Map<String, Object> userInputs = Map.of("file1", "gs://bucket3/path/to/file.vcf.gz");

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    Set<String> uniqueBucketsResult =
        pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
            inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun);

    assertEquals(0, uniqueBucketsResult.size());
  }

  @Test
  void extractUniqueBucketsFromManifestsMissingOptionalManifestOk() {
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true),
            createTestPipelineInputDefWithName(
                "manifest2", "manifest_2", PipelineVariableTypesEnum.MANIFEST, false, true));

    String manifestFile1 = "gs://bucket1/path/to/manifest1.tsv";

    String file1 = "gs://bucket4/path/to/file1.vcf.gz";
    String file2 = "gs://bucket5/path/to/file2.vcf.gz";

    String manifestFile1Contents = "sample1\t%s%nsample2\t%s%n".formatted(file1, file2);

    // expected buckets are the buckets from the files, not from the manifests
    Set<String> expectedBucketSet = Set.of("bucket4", "bucket5");

    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile1Contents));

    Map<String, Object> userInputs = Map.of("manifest1", manifestFile1);

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    Set<String> uniqueBucketsResult =
        pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
            inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun);

    assertEquals(expectedBucketSet.size(), uniqueBucketsResult.size());
    assertEquals(expectedBucketSet, uniqueBucketsResult);
  }

  @Test
  void extractUniqueBucketsFromManifestsNoGcsFiles() {
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true));

    String manifestFile1 = "gs://bucket1/path/to/manifest1.tsv";

    String file1 = "file1.vcf.gz";
    String file2 = "file2.vcf.gz";

    String manifestFile1Contents = "sample1\t%s%nsample2\t%s%n".formatted(file1, file2);

    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile1Contents));

    Map<String, Object> userInputs = Map.of("manifest1", manifestFile1);

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    assertThrows(
        ValidationException.class,
        () ->
            pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
                inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun));
  }

  @Test
  void extractUniqueBucketsFromManifestsEmptyLineInMiddle() {
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true));

    String manifestFile1 = "gs://bucket1/path/to/manifest1.tsv";

    String file1 = "file1.vcf.gz";
    String file2 = "file2.vcf.gz";

    String manifestFile1Contents = "sample1\t%s%n%nsample2\t%s%n".formatted(file1, file2);

    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile1Contents));

    Map<String, Object> userInputs = Map.of("manifest1", manifestFile1);

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    assertThrows(
        ValidationException.class,
        () ->
            pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
                inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun));
  }

  @Test
  void extractUniqueBucketsFromManifestsLocal() {
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true),
            createTestPipelineInputDefWithName(
                "manifest2", "manifest_2", PipelineVariableTypesEnum.MANIFEST, false, true),
            createTestPipelineInputDefWithName(
                "file1", "file_1", PipelineVariableTypesEnum.FILE, false, true));

    String manifestFile1 = "path/to/manifest1.tsv";
    String manifestFile2 = "path/to/manifest2.tsv";
    String manifestFile1FullPath =
        "gs://%s/%s"
            .formatted(
                TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
                constructDestinationBlobNameForUserInputFile(TEST_JOB_ID, manifestFile1));
    String manifestFile2FullPath =
        "gs://%s/%s"
            .formatted(
                TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
                constructDestinationBlobNameForUserInputFile(TEST_JOB_ID, manifestFile2));

    String file1 = "gs://bucket4/path/to/file1.vcf.gz";
    String file2 = "gs://bucket5/path/to/file2.vcf.gz";
    String file3 = "gs://bucket6/path/to/file3.vcf.gz";
    String file4 = "gs://bucket6/path/to/file4.vcf.gz";

    String manifestFile1Contents = "sample1\t%s%nsample2\t%s%n".formatted(file1, file2);
    String manifestFile2Contents = "sample3\t%s%nsample4\t%s%n".formatted(file3, file4);

    // expected buckets are the buckets from the files, not from the manifests
    Set<String> expectedBucketSet = Set.of("bucket4", "bucket5", "bucket6");

    // local inputs will have been copied into control workspace bucket
    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1FullPath)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile1Contents));
    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile2FullPath)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile2Contents));

    Map<String, Object> userInputs =
        Map.of(
            "manifest1", manifestFile1,
            "manifest2", manifestFile2,
            "file1", "gs://bucket3/path/to/file.vcf.gz");

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    Set<String> uniqueBucketsResult =
        pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
            inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun);

    assertEquals(expectedBucketSet.size(), uniqueBucketsResult.size());
    assertEquals(expectedBucketSet, uniqueBucketsResult);
  }

  @Test
  void extractUniqueBucketsFromManifestsThrowsException() {
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true));

    String manifestFile1 = "gs://bucket1/path/to/manifest1.tsv";

    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1)))
        .thenThrow(new RuntimeException("could not read manifest file"));

    Map<String, Object> userInputs = Map.of("manifest1", manifestFile1);

    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    // we rethrow any exception as an InternalServerErrorException
    assertThrows(
        InternalServerErrorException.class,
        () ->
            pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
                inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun));
  }

  @Test
  void extractUniqueBucketsFromManifestsMismatchedItemsPerLineThrowsValidationException() {
    List<PipelineInputDefinition> inputDefinitions =
        List.of(
            createTestPipelineInputDefWithName(
                "manifest1", "manifest_1", PipelineVariableTypesEnum.MANIFEST, false, true));

    String manifestFile1 = "gs://bucket1/path/to/manifest1.tsv";
    String file1 = "gs://bucket4/path/to/file1.vcf.gz";
    String file2 = "gs://bucket5/path/to/file2.vcf.gz";

    // First row has 2 items, second row has 3 items.
    String manifestFile1Contents =
        "sample1\t%s%nsample2\t%s\textra-column%n".formatted(file1, file2);

    when(mockGcsService.getBufferedReaderForGcsTextFile(new GcsFile(manifestFile1)))
        .thenReturn(getBufferedReaderForStringTesting(manifestFile1Contents));

    Map<String, Object> userInputs = Map.of("manifest1", manifestFile1);
    PipelineRun pipelineRun = createAndSavePipelineRunWithInputs(userInputs);

    ValidationException e =
        assertThrows(
            ValidationException.class,
            () ->
                pipelineInputsOutputsService.extractUniqueBucketsFromManifests(
                    inputDefinitions, TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME, pipelineRun));
    assertTrue(
        e.getMessage()
            .contains(
                "Manifest file manifest1.tsv has inconsistent number of items at line 2. Expected 2 items, found 3."));
  }

  // helper method to create and save a pipeline run with the given user inputs
  PipelineRun createAndSavePipelineRunWithInputs(Map<String, Object> userInputs) {
    PipelineRun pipelineRun = TestUtils.createNewPipelineRunWithJobId(TEST_JOB_ID);
    pipelineRunsRepository.save(pipelineRun);

    // store user inputs in the db
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setInputs(pipelineInputsOutputsService.mapToString(userInputs));
    pipelineInput.setPipelineRunId(pipelineRun.getId());
    pipelineInputsRepository.save(pipelineInput);

    return pipelineRun;
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
                TestUtils.createTestPipelineInputDefWithName(
                    "inputName",
                    "input_name",
                    PipelineVariableTypesEnum.STRING,
                    false,
                    true,
                    false,
                    null,
                    null,
                    null),
                TestUtils.createTestPipelineInputDefWithName(
                    "inputNameOptionalNotSpecified",
                    "input_name_optional_not_specified",
                    PipelineVariableTypesEnum.STRING,
                    false,
                    true,
                    false,
                    "default optional value not specified by user",
                    null,
                    null),
                TestUtils.createTestPipelineInputDefWithName(
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
                TestUtils.createTestPipelineInputDefWithName(
                    "inputName",
                    "input_name",
                    PipelineVariableTypesEnum.STRING,
                    true,
                    true,
                    false,
                    null,
                    null,
                    null),
                TestUtils.createTestPipelineInputDefWithName(
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
                TestUtils.createTestPipelineInputDefWithName(
                    "inputNameUserProvided",
                    "input_name_user_provided",
                    PipelineVariableTypesEnum.STRING,
                    true,
                    true,
                    false,
                    null,
                    null,
                    null),
                TestUtils.createTestPipelineInputDefWithName(
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
            TestUtils.createTestPipelineInputDefWithName(
                "userInputName",
                "user_input_name",
                PipelineVariableTypesEnum.FILE,
                true,
                true,
                false,
                null,
                null,
                null),
            TestUtils.createTestPipelineInputDefWithName(
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

  @Test
  void savePipelineOutputsWithFileSizesSuccess() {
    PipelineRun newPipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
    pipelineRunsRepository.save(newPipelineRun);

    pipelineInputsOutputsService.savePipelineOutputs(
        newPipelineRun.getId(),
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE,
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE_SIZE);

    List<PipelineOutput> savedOutputs =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(newPipelineRun.getId());
    assertEquals(2, savedOutputs.size());

    PipelineOutput fileOutput =
        savedOutputs.stream()
            .filter(o -> o.getOutputName().equals("testFileOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals(
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.get("testFileOutputKey"),
        fileOutput.getOutputValue());
    assertEquals(256L, fileOutput.getFileSizeBytes());

    PipelineOutput stringOutput =
        savedOutputs.stream()
            .filter(o -> o.getOutputName().equals("testStringOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals(
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.get("testStringOutputKey"),
        stringOutput.getOutputValue());
    assertNull(stringOutput.getFileSizeBytes());
  }

  @Test
  void savePipelineOutputsWithoutFileSizesSuccess() {
    PipelineRun newPipelineRun = createNewPipelineRunWithJobId(UUID.randomUUID());
    pipelineRunsRepository.save(newPipelineRun);

    pipelineInputsOutputsService.savePipelineOutputs(
        newPipelineRun.getId(), TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE, Collections.emptyMap());

    List<PipelineOutput> savedOutputs =
        pipelineOutputsRepository.findPipelineOutputsByPipelineRunId(newPipelineRun.getId());
    assertEquals(2, savedOutputs.size());

    PipelineOutput fileOutput =
        savedOutputs.stream()
            .filter(o -> o.getOutputName().equals("testFileOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals(
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.get("testFileOutputKey"),
        fileOutput.getOutputValue());
    assertNull(fileOutput.getFileSizeBytes());

    PipelineOutput stringOutput =
        savedOutputs.stream()
            .filter(o -> o.getOutputName().equals("testStringOutputKey"))
            .findFirst()
            .orElseThrow();
    assertEquals(
        TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.get("testStringOutputKey"),
        stringOutput.getOutputValue());
    assertNull(stringOutput.getFileSizeBytes());
  }

  @Test
  void getPipelineOutputsFileSizeSuccess() {
    Pipeline testPipeline = updateTestPipeline1WithTestValues();
    Long pipelineId = testPipeline.getId();

    PipelineOutputDefinition fileOutput1 =
        new PipelineOutputDefinition(
            pipelineId,
            "fileOutput1",
            "file_output_1",
            null,
            null,
            PipelineVariableTypesEnum.FILE,
            true);
    PipelineOutputDefinition fileOutput2 =
        new PipelineOutputDefinition(
            pipelineId,
            "fileOutput2",
            "file_output_2",
            null,
            null,
            PipelineVariableTypesEnum.FILE,
            true);
    PipelineOutputDefinition stringOutput =
        new PipelineOutputDefinition(
            pipelineId,
            "stringOutput",
            "string_output",
            null,
            null,
            PipelineVariableTypesEnum.STRING,
            true);

    testPipeline.setPipelineOutputDefinitions(List.of(fileOutput1, fileOutput2, stringOutput));

    String filePath1 = "gs://bucket/path/file1.vcf.gz";
    String filePath2 = "gs://bucket/path/file2.vcf.gz";
    Long fileSize1 = 12345L;
    Long fileSize2 = 6789L;

    Map<String, String> outputsMap =
        Map.of(
            "fileOutput1", filePath1,
            "fileOutput2", filePath2,
            "stringOutput", "IAmGroot");

    when(mockGcsService.getFileSizeInBytes(filePath1)).thenReturn(fileSize1);
    when(mockGcsService.getFileSizeInBytes(filePath2)).thenReturn(fileSize2);

    Map<String, Long> outputFileSizes =
        pipelineInputsOutputsService.getPipelineOutputsFileSize(testPipeline, outputsMap);

    assertEquals(2, outputFileSizes.size());
    assertEquals(fileSize1, outputFileSizes.get("fileOutput1"));
    assertEquals(fileSize2, outputFileSizes.get("fileOutput2"));
  }

  @Test
  void getPipelineOutputsFileSizeMissingOutputThrowsException() {
    Pipeline testPipeline = updateTestPipeline1WithTestValues();
    testPipeline.setPipelineOutputDefinitions(TestUtils.TEST_PIPELINE_OUTPUT_DEFINITIONS_WITH_FILE);

    String fileOutputKey = "testFileOutputKey";
    Map<String, String> outputsMap = Map.of("testStringOutputKey", "IAmGroot");

    InternalServerErrorException exception =
        assertThrows(
            InternalServerErrorException.class,
            () ->
                pipelineInputsOutputsService.getPipelineOutputsFileSize(testPipeline, outputsMap));

    assertEquals(
        "File output %s is missing from outputs map".formatted(fileOutputKey),
        exception.getMessage());
  }

  // test helper methods
  private static PipelineInputDefinition createTestPipelineInputDef(
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided,
      boolean isCustomValue,
      String defaultValue) {
    return createTestPipelineInputDefWithNameAndPipelineId(
        3L,
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
    return createTestPipelineInputDefWithNameAndPipelineId(
        3L,
        inputName,
        inputWdlVariableName,
        type,
        isRequired,
        isUserProvided,
        false,
        null,
        null,
        null);
  }

  private static PipelineInputDefinition createTestPipelineInputDefWithNameAndPipelineId(
      Long pipelineId,
      String inputName,
      String inputWdlVariableName,
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided) {
    return createTestPipelineInputDefWithNameAndPipelineId(
        pipelineId,
        inputName,
        inputWdlVariableName,
        type,
        isRequired,
        isUserProvided,
        false,
        null,
        null,
        null);
  }

  private static PipelineInputDefinition createTestPipelineInputDefWithNameAndPipelineId(
      Long pipelineId,
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
        switch (type) {
          case FILE, FILE_ARRAY -> ".vcf.gz";
          case MANIFEST -> ".tsv";
          default -> null;
        };
    return new PipelineInputDefinition(
        pipelineId,
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

  @Test
  void deleteOutputSourcesFilesSuccess() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);

    Map<String, Object> outputsMap = new HashMap<>();
    outputsMap.put("outputFile1", "gs://source-bucket/path/to/file1.vcf.gz");
    outputsMap.put("outputFile2", "gs://source-bucket/path/to/file2.vcf.gz");

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    ArgumentCaptor<GcsFile> fileCaptor = ArgumentCaptor.forClass(GcsFile.class);

    doNothing().when(mockGcsService).deleteObject(any(GcsFile.class));

    pipelineInputsOutputsService.deleteOutputSourcesFiles(testPipelineRun);

    // Verify deleteObject was called twice (once for each file)
    verify(mockGcsService, times(2)).deleteObject(fileCaptor.capture());

    // Verify the files that were deleted
    List<GcsFile> deletedFiles = fileCaptor.getAllValues();
    assertEquals(2, deletedFiles.size());
    assertTrue(
        deletedFiles.stream()
            .anyMatch(
                file -> file.getFullPath().equals("gs://source-bucket/path/to/file1.vcf.gz")));
    assertTrue(
        deletedFiles.stream()
            .anyMatch(
                file -> file.getFullPath().equals("gs://source-bucket/path/to/file2.vcf.gz")));
  }

  @Test
  void deleteOutputSourcesFilesWithDeletionFailureDoesNotThrow() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);

    Map<String, Object> outputsMap = new HashMap<>();
    outputsMap.put("outputFile", "gs://source-bucket/path/to/file.vcf.gz");

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    // Mock GCS service to throw an exception
    doThrow(new RuntimeException("Failed to delete file"))
        .when(mockGcsService)
        .deleteObject(any(GcsFile.class));

    // Should not throw an exception even when deletion fails (just logs error)
    assertDoesNotThrow(
        () -> pipelineInputsOutputsService.deleteOutputSourcesFiles(testPipelineRun));

    // Verify that deleteObject was still called
    verify(mockGcsService).deleteObject(any(GcsFile.class));
  }

  @Test
  void deleteOutputSourcesFilesWithEmptyOutputsMap() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);

    Map<String, Object> outputsMap = new HashMap<>();

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    // Should not throw an exception even with empty outputs
    assertDoesNotThrow(
        () -> pipelineInputsOutputsService.deleteOutputSourcesFiles(testPipelineRun));

    // Verify that deleteObject was never called
    verify(mockGcsService, never()).deleteObject(any(GcsFile.class));
  }

  @Test
  void deleteOutputSourcesFilesWithMultipleFilesAndPartialFailure() {
    Pipeline testPipeline = createTestPipelineWithId();
    UUID jobId = UUID.randomUUID();
    PipelineRun testPipelineRun = createTestPipelineRun(testPipeline, jobId);

    Map<String, Object> outputsMap = new HashMap<>();
    outputsMap.put("outputFile1", "gs://source-bucket/path/to/file1.vcf.gz");
    outputsMap.put("outputFile2", "gs://source-bucket/path/to/file2.vcf.gz");
    outputsMap.put("outputFile3", "gs://source-bucket/path/to/file3.vcf.gz");

    PipelineOutput pipelineOutput = new PipelineOutput();
    pipelineOutput.setJobId(testPipelineRun.getId());
    pipelineOutput.setOutputs(pipelineInputsOutputsService.mapToString(outputsMap));
    pipelineOutputsRepository.save(pipelineOutput);

    // Mock GCS service to fail on second file only
    doNothing()
        .doThrow(new RuntimeException("Failed to delete file2"))
        .doNothing()
        .when(mockGcsService)
        .deleteObject(any(GcsFile.class));

    // Should not throw an exception even when one deletion fails
    assertDoesNotThrow(
        () -> pipelineInputsOutputsService.deleteOutputSourcesFiles(testPipelineRun));

    // Verify that deleteObject was called three times (attempts all deletions despite failure)
    verify(mockGcsService, times(3)).deleteObject(any(GcsFile.class));
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

  private static List<PipelineOutput> getPipelineOutputsForPipelineRun(
      PipelineRun pipelineRun, Boolean storeFileSize) {
    return TestUtils.TEST_PIPELINE_OUTPUTS_WITH_FILE.entrySet().stream()
        .map(
            entry -> {
              PipelineOutput pipelineOutput = new PipelineOutput();
              pipelineOutput.setPipelineRunId(pipelineRun.getId());
              pipelineOutput.setOutputName(entry.getKey());
              pipelineOutput.setOutputValue(entry.getValue());
              // Only set file size for file outputs
              if (storeFileSize && entry.getKey().equals("testFileOutputKey")) {
                pipelineOutput.setFileSizeBytes(256L);
              }
              return pipelineOutput;
            })
        .toList();
  }
}
