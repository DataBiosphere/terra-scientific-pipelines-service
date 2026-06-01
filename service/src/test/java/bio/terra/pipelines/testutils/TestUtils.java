package bio.terra.pipelines.testutils;

import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import bio.terra.pipelines.model.PipelineQuota;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.rawls.model.MethodConfiguration;
import bio.terra.rawls.model.MethodRepoMethod;
import java.io.BufferedReader;
import java.io.StringReader;
import java.math.BigDecimal;
import java.util.*;
import org.mockito.ArgumentMatcher;

/** A collection of utilities and constants useful for tests. */
public class TestUtils {
  // Pipelines test constants
  public static final PipelinesEnum TEST_PIPELINE_1_IMPUTATION_ENUM =
      PipelinesEnum.ARRAY_IMPUTATION;

  public static final int TEST_PIPELINE_VERSION_1 = 1;
  public static final String TEST_PIPELINE_KEY_1 = "array_imputation_v1";
  public static final boolean TEST_PIPELINE_HIDDEN_1 = false;
  public static final String TEST_PIPELINE_DISPLAY_NAME_1 =
      "Test Pipeline Name"; // this matches the job pre-populated in the db for tests
  public static final String TEST_PIPELINE_DESCRIPTION_1 = "Test Pipeline Description";
  public static final String TEST_PIPELINE_TYPE_1 = "imputation1";
  public static final String TEST_TOOL_NAME_1 = "methodName1";
  public static final String TEST_TOOL_NAME_WITH_PIPELINE_VERSION_1 = "methodName1_v1";
  public static final String TEST_DATA_TABLE_ENTITY_NAME_1 = "array_imputation_v1";
  public static final String TEST_TOOL_VERSION_1 = "0.1.12";
  public static final int TEST_PIPELINE_VERSION_2 = 1;
  public static final boolean TEST_PIPELINE_HIDDEN_2 = false;
  public static final String TEST_PIPELINE_DISPLAY_NAME_2 = "Test Pipeline Name Two";
  public static final String TEST_PIPELINE_DESCRIPTION_2 = "Test Pipeline Description Two";
  public static final String TEST_PIPELINE_TYPE_2 = "imputation2";
  public static final String TEST_TOOL_NAME_2 = "methodName2";
  public static final String TEST_TOOL_VERSION_2 = "1.1.12";
  public static final int TEST_LOW_PASS_IMPUTATION_PIPELINE_VERSION = 1;
  public static final UUID CONTROL_WORKSPACE_ID =
      UUID.fromString("fafafafa-fafa-fafa-fafa-fafafafafafa");
  public static final String CONTROL_WORKSPACE_BILLING_PROJECT = "testTerraProject";
  public static final String CONTROL_WORKSPACE_NAME = "testTerraWorkspaceName";
  public static final String CONTROL_WORKSPACE_CONTAINER_NAME =
      "fc-secure-%s".formatted(CONTROL_WORKSPACE_ID);
  public static final String GCP_STORAGE_PROTOCOL = "gs://";
  public static final String CONTROL_WORKSPACE_GOOGLE_PROJECT = "testGoogleProject";
  public static final List<PipelineOutputDefinition> TEST_PIPELINE_OUTPUT_DEFINITIONS_WITH_FILE =
      new ArrayList<>(
          List.of(
              PipelineOutputDefinition.builder()
                  .name("testFileOutputKey")
                  .wdlVariableName("test_file_output_key")
                  .displayName("Test File Output Display Name")
                  .description("test output file description")
                  .type(PipelineVariableTypesEnum.FILE)
                  .isRequired(true)
                  .build(),
              PipelineOutputDefinition.builder()
                  .name("testStringOutputKey")
                  .wdlVariableName("test_string_output_key")
                  .displayName("Test String Output Display Name")
                  .description("test output string description")
                  .type(PipelineVariableTypesEnum.STRING)
                  .isRequired(true)
                  .build()));

  public static final Map<String, String> TEST_PIPELINE_OUTPUTS_WITH_FILE =
      new HashMap(
          Map.of(
              "testFileOutputKey",
              "gs://fc-secure-%s/testFileOutputValue".formatted(CONTROL_WORKSPACE_ID),
              "testStringOutputKey",
              "testStringOutputValue"));
  public static final Map<String, Long> TEST_PIPELINE_OUTPUTS_WITH_FILE_SIZE =
      new HashMap(Map.of("testFileOutputKey", 256L));
  public static final Map<String, Object> TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED_V3 =
      new HashMap<>(
          Map.of(
              "testFileOutputKey",
              Map.of("value", "testFileOutputValue"),
              "testStringOutputKey",
              Map.of("value", "testStringOutputValue")));

  public static final Map<String, String> TEST_PIPELINE_OUTPUTS_WITH_FILE_FORMATTED =
      new HashMap(
          Map.of(
              "testFileOutputKey",
              "testFileOutputValue",
              "testStringOutputKey",
              "testStringOutputValue"));

  public static final List<PipelineInputDefinition> TEST_PIPELINE_INPUTS_DEFINITION_LIST =
      new ArrayList<>(
          List.of(
              PipelineInputDefinition.builder()
                  .name("testRequiredStringInput")
                  .wdlVariableName("test_required_string_input")
                  .displayName("test required string input's display name")
                  .description("the description of this input")
                  .type(PipelineVariableTypesEnum.STRING)
                  .fileSuffix(null)
                  .isRequired(true)
                  .userProvided(true)
                  .defaultValue(null)
                  .minValue(null)
                  .maxValue(null)
                  .build(),
              PipelineInputDefinition.builder()
                  .name("testOptionalStringInput")
                  .wdlVariableName("test_optional_string_input")
                  .displayName("test optional string input's display name")
                  .description("the description of this optional input")
                  .type(PipelineVariableTypesEnum.STRING)
                  .fileSuffix(null)
                  .isRequired(false)
                  .userProvided(true)
                  .defaultValue("testDefaultValue")
                  .minValue(null)
                  .maxValue(null)
                  .build(),
              PipelineInputDefinition.builder()
                  .name("testRequiredIntInput")
                  .wdlVariableName("test_required_int_input")
                  .displayName("test required int input's display name")
                  .description("the description of this int input")
                  .type(PipelineVariableTypesEnum.INTEGER)
                  .fileSuffix(null)
                  .isRequired(true)
                  .userProvided(true)
                  .defaultValue(null)
                  .minValue(0.0)
                  .maxValue(1.0)
                  .build(),
              PipelineInputDefinition.builder()
                  .name("testOptionalIntInput")
                  .wdlVariableName("test_optional_int_input")
                  .displayName("test optional int input's display name")
                  .description("the description of this optional int input")
                  .type(PipelineVariableTypesEnum.INTEGER)
                  .fileSuffix(null)
                  .isRequired(false)
                  .userProvided(true)
                  .defaultValue("42")
                  .minValue(0.0)
                  .maxValue(100.0)
                  .build(),
              PipelineInputDefinition.builder()
                  .name("testServiceProvidedInput")
                  .wdlVariableName("test_service_provided_input")
                  .displayName(null)
                  .description(null)
                  .type(PipelineVariableTypesEnum.STRING)
                  .fileSuffix(null)
                  .isRequired(true)
                  .userProvided(false)
                  .defaultValue("testServiceProvidedDefaultValue")
                  .minValue(null)
                  .maxValue(null)
                  .build(),
              PipelineInputDefinition.builder()
                  .name("testRequiredVcfInput")
                  .wdlVariableName("test_required_vcf_input")
                  .displayName(null)
                  .description(null)
                  .type(PipelineVariableTypesEnum.FILE)
                  .fileSuffix(".vcf.gz")
                  .isRequired(true)
                  .userProvided(true)
                  .defaultValue(null)
                  .minValue(null)
                  .maxValue(null)
                  .build()));

  public static final List<PipelineOutputDefinition> TEST_PIPELINE_OUTPUTS_DEFINITION_LIST =
      new ArrayList<>(
          List.of(
              PipelineOutputDefinition.builder()
                  .name("outputString")
                  .wdlVariableName("output_string")
                  .displayName("output string display name")
                  .description("description")
                  .type(PipelineVariableTypesEnum.STRING)
                  .isRequired(true)
                  .build(),
              PipelineOutputDefinition.builder()
                  .name("outputInteger")
                  .wdlVariableName("output_integer")
                  .displayName("output integer display name")
                  .description("description")
                  .type(PipelineVariableTypesEnum.INTEGER)
                  .isRequired(true)
                  .build(),
              PipelineOutputDefinition.builder()
                  .name("outputBooleanOptional")
                  .wdlVariableName("output_boolean_optional")
                  .displayName(null)
                  .description(null)
                  .type(PipelineVariableTypesEnum.BOOLEAN)
                  .isRequired(false)
                  .build()));

  public static final Map<String, Object> TEST_PIPELINE_OUTPUTS_FROM_ENTITY =
      new HashMap<>(
          Map.of(
              "output_string",
              "test",
              "output_integer",
              123,
              "output_boolean_optional",
              false)); // matches TestUtils.TOOL_CONFIG_GENERIC output definitions
  public static final Map<String, Object> TEST_PIPELINE_OUTPUTS_PROCESSED =
      new HashMap<>(
          Map.of(
              "outputString",
              "test",
              "outputInteger",
              123,
              "outputBooleanOptional",
              false)); // matches TestUtils.TOOL_CONFIG_GENERIC output definitions
  public static final Pipeline TEST_ARRAY_IMPUTATION_PIPELINE_1 =
      Pipeline.builder()
          .name(PipelinesEnum.ARRAY_IMPUTATION)
          .version(TEST_PIPELINE_VERSION_1)
          .pipelineKey(buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, TEST_PIPELINE_VERSION_1))
          .hidden(TEST_PIPELINE_HIDDEN_1)
          .displayName(TEST_PIPELINE_DISPLAY_NAME_1)
          .description(TEST_PIPELINE_DESCRIPTION_1)
          .pipelineType(TEST_PIPELINE_TYPE_1)
          .toolName(TEST_TOOL_NAME_1)
          .toolVersion(TEST_TOOL_VERSION_1)
          .workspaceBillingProject(CONTROL_WORKSPACE_BILLING_PROJECT)
          .workspaceName(CONTROL_WORKSPACE_NAME)
          .workspaceStorageContainerName(CONTROL_WORKSPACE_CONTAINER_NAME)
          .workspaceGoogleProject(CONTROL_WORKSPACE_GOOGLE_PROJECT)
          .inputDefinitions(TEST_PIPELINE_INPUTS_DEFINITION_LIST)
          .outputDefinitions(TEST_PIPELINE_OUTPUTS_DEFINITION_LIST)
          .build();

  public static final Pipeline TEST_ARRAY_IMPUTATION_PIPELINE_2 =
      Pipeline.builder()
          .name(PipelinesEnum.ARRAY_IMPUTATION)
          .version(TEST_PIPELINE_VERSION_2)
          .pipelineKey(buildPipelineKey(PipelinesEnum.ARRAY_IMPUTATION, TEST_PIPELINE_VERSION_2))
          .hidden(TEST_PIPELINE_HIDDEN_2)
          .displayName(TEST_PIPELINE_DISPLAY_NAME_2)
          .description(TEST_PIPELINE_DESCRIPTION_2)
          .pipelineType(TEST_PIPELINE_TYPE_2)
          .toolName(TEST_TOOL_NAME_2)
          .toolVersion(TEST_TOOL_VERSION_2)
          .workspaceBillingProject(CONTROL_WORKSPACE_BILLING_PROJECT)
          .workspaceName(CONTROL_WORKSPACE_NAME)
          .workspaceStorageContainerName(CONTROL_WORKSPACE_CONTAINER_NAME)
          .workspaceGoogleProject(CONTROL_WORKSPACE_GOOGLE_PROJECT)
          .inputDefinitions(TEST_PIPELINE_INPUTS_DEFINITION_LIST)
          .outputDefinitions(TEST_PIPELINE_OUTPUTS_DEFINITION_LIST)
          .build();

  public static final Pipeline TEST_LOW_PASS_IMPUTATION_PIPELINE =
      Pipeline.builder()
          .name(PipelinesEnum.LOW_PASS_IMPUTATION)
          .version(TEST_LOW_PASS_IMPUTATION_PIPELINE_VERSION)
          .pipelineKey(
              buildPipelineKey(
                  PipelinesEnum.LOW_PASS_IMPUTATION, TEST_LOW_PASS_IMPUTATION_PIPELINE_VERSION))
          .hidden(TEST_PIPELINE_HIDDEN_1)
          .displayName(TEST_PIPELINE_DISPLAY_NAME_1)
          .description(TEST_PIPELINE_DESCRIPTION_1)
          .pipelineType(TEST_PIPELINE_TYPE_1)
          .toolName(TEST_TOOL_NAME_1)
          .toolVersion(TEST_TOOL_VERSION_1)
          .workspaceBillingProject(CONTROL_WORKSPACE_BILLING_PROJECT)
          .workspaceName(CONTROL_WORKSPACE_NAME)
          .workspaceStorageContainerName(CONTROL_WORKSPACE_CONTAINER_NAME)
          .workspaceGoogleProject(CONTROL_WORKSPACE_GOOGLE_PROJECT)
          .inputDefinitions(TEST_PIPELINE_INPUTS_DEFINITION_LIST)
          .outputDefinitions(TEST_PIPELINE_OUTPUTS_DEFINITION_LIST)
          .build();

  public static final PipelineQuota TEST_PIPELINE_QUOTA_1 =
      PipelineQuota.builder()
          .defaultQuota(100)
          .minQuotaConsumed(10)
          .quotaUnits(QuotaUnitsEnum.SAMPLES)
          .build();

  public static final String TEST_USER_1_ID =
      "testUser"; // this matches the job pre-populated in the db for tests
  public static final String TEST_USER_1_EMAIL = "testUser@test.com";
  public static final BearerToken TEST_USER_1_BEARER_TOKEN = new BearerToken("fake-token");
  public static final SamUser TEST_SAM_USER_1 =
      new SamUser(TEST_USER_1_EMAIL, TEST_USER_1_ID, TEST_USER_1_BEARER_TOKEN);
  public static final String TEST_USER_2_ID = "testUser2";

  public static final UUID TEST_NEW_UUID = UUID.fromString("deadbeef-dead-beef-aaaa-beefdeadbeef");

  public static final UUID TEST_NEW_UUID_2 =
      UUID.fromString("deadbeef-dead-beef-bbbb-beefdeadbeef");

  public static final Map<String, Object> TEST_PIPELINE_INPUTS =
      new HashMap<>(Map.of("first_key", "first_value"));

  public static final Map<String, Object> TEST_PIPELINE_INPUTS_ARRAY_IMPUTATION =
      new HashMap<>(
          Map.of("multiSampleVcf", "fake/file.vcf.gz", "outputBasename", "fake_basename"));

  public static final String TEST_DOMAIN = "some-teaspoons-domain.com";
  public static final String TEST_USER_PROVIDED_DESCRIPTION =
      "user-provided description of a pipeline run";

  // tool config stuff
  public static final String TOOL_CONFIG_KEY = "tool_config_key";
  public static final String TOOL_OUTPUTS_KEY = "tool_outputs_key";
  public static final String SUBMISSION_ID_KEY = "submission_id_key";
  public static final ToolConfig TOOL_CONFIG_GENERIC =
      new ToolConfig(
          TestUtils.TEST_TOOL_NAME_1,
          TestUtils.TEST_TOOL_VERSION_1,
          TestUtils.TEST_TOOL_NAME_WITH_PIPELINE_VERSION_1,
          TestUtils.TEST_DATA_TABLE_ENTITY_NAME_1,
          TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST,
          TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST,
          true,
          "gs://path/to/monitoring/script.sh",
          false,
          BigDecimal.valueOf(2.0),
          1L);

  public static final MethodConfiguration VALID_METHOD_CONFIGURATION =
      new MethodConfiguration()
          .name("name")
          .inputs(Map.of("workflowName.first_input", "this.first_input"))
          .outputs(Map.of("workflowName.first_output", "this.first_output"))
          .methodRepoMethod(
              new MethodRepoMethod()
                  .methodName("methodName")
                  .methodNamespace("namespace")
                  .methodVersion("1.2.3")
                  .methodUri("this/is/a/uri/with/a/version/1.2.3"))
          .rootEntityType("imputation_beagle");

  public static PipelineRun createNewPipelineRunWithJobId(UUID jobId) {
    return createNewPipelineRunWithJobIdAndUser(jobId, TEST_USER_1_ID);
  }

  public static PipelineRun createNewPipelineRunWithJobIdAndUser(UUID jobId, String userId) {
    return new PipelineRun(
        jobId,
        userId,
        TEST_PIPELINE_KEY_1,
        TEST_TOOL_VERSION_1,
        CONTROL_WORKSPACE_BILLING_PROJECT,
        CONTROL_WORKSPACE_NAME,
        CONTROL_WORKSPACE_CONTAINER_NAME,
        CONTROL_WORKSPACE_GOOGLE_PROJECT,
        CommonPipelineRunStatusEnum.PREPARING,
        "test description");
  }

  public static String buildTestResultUrl(String jobId, int resultApiVersion) {
    return "https://%s/api/pipelineruns/v%s/result/%s"
        .formatted(TEST_DOMAIN, resultApiVersion, jobId);
  }

  /**
   * Creates a test pipeline with the specified parameters. Uses default values for common fields
   * like workspace details.
   */
  public static Pipeline createTestPipeline(
      PipelinesEnum name, int version, boolean hidden, String displayName, String toolVersion) {
    return Pipeline.builder()
        .name(name)
        .version(version)
        .pipelineKey(buildPipelineKey(name, version))
        .hidden(hidden)
        .displayName(displayName)
        .description("description")
        .pipelineType("pipelineType")
        .toolName("toolName")
        .toolVersion(toolVersion)
        .workspaceBillingProject(CONTROL_WORKSPACE_BILLING_PROJECT)
        .workspaceName(CONTROL_WORKSPACE_NAME)
        .workspaceStorageContainerName(CONTROL_WORKSPACE_CONTAINER_NAME)
        .workspaceGoogleProject(CONTROL_WORKSPACE_GOOGLE_PROJECT)
        .build();
  }

  public static String buildPipelineKey(PipelinesEnum pipelineName, Integer pipelineVersion) {
    return "%s_v%s".formatted(pipelineName.getValue(), pipelineVersion);
  }

  /** Helper method to create a BufferedReader from a string for testing purposes. */
  public static BufferedReader getBufferedReaderForStringTesting(String fileContentsString) {
    return new BufferedReader(new StringReader(fileContentsString));
  }

  public static PipelineInputDefinition createTestPipelineInputDef(
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided,
      boolean isCustomValue,
      String defaultValue) {
    return createTestPipelineInputDefWithName(
        "inputName", "input_name", type, isRequired, isUserProvided, defaultValue, null, null);
  }

  public static PipelineInputDefinition createTestPipelineInputDefWithName(
      String inputName,
      String inputWdlVariableName,
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided) {
    return createTestPipelineInputDefWithName(
        inputName, inputWdlVariableName, type, isRequired, isUserProvided, null, null, null);
  }

  public static PipelineInputDefinition createTestPipelineInputDefWithName(
      String inputName,
      String inputWdlVariableName,
      PipelineVariableTypesEnum type,
      boolean isRequired,
      boolean isUserProvided,
      String defaultValue,
      Double minValue,
      Double maxValue) {
    String fileSuffix =
        switch (type) {
          case FILE, FILE_ARRAY -> ".vcf.gz";
          case MANIFEST -> ".tsv";
          default -> null;
        };
    return PipelineInputDefinition.builder()
        .name(inputName)
        .wdlVariableName(inputWdlVariableName)
        .displayName(null)
        .description(null)
        .type(type)
        .fileSuffix(fileSuffix)
        .isRequired(isRequired)
        .userProvided(isUserProvided)
        .defaultValue(defaultValue)
        .minValue(minValue)
        .maxValue(maxValue)
        .build();
  }

  /**
   * Helper method to assess "equality" of lists of input definitions when the instances of the list
   * or its elements are different but the content is the same. This is useful for matching
   * arguments in Mockito.
   *
   * @param expectedInputDefinitions
   * @return
   */
  public static ArgumentMatcher<List<PipelineInputDefinition>> matchesExpectedInputDefinitions(
      List<PipelineInputDefinition> expectedInputDefinitions) {
    List<String> expectedNames =
        expectedInputDefinitions.stream().map(PipelineInputDefinition::getName).toList();

    return defs ->
        defs != null
            && defs.size() == expectedNames.size()
            && defs.stream().map(PipelineInputDefinition::getName).toList().equals(expectedNames);
  }

  /**
   * Helper method to assess "equality" of lists of output definitions when the instances of the
   * list or its elements are different but the content is the same. This is useful for matching
   * arguments in Mockito.
   *
   * @param expectedOutputDefinitions
   * @return
   */
  public static ArgumentMatcher<List<PipelineOutputDefinition>> matchesExpectedOutputDefinitions(
      List<PipelineOutputDefinition> expectedOutputDefinitions) {
    List<String> expectedNames =
        expectedOutputDefinitions.stream().map(PipelineOutputDefinition::getName).toList();

    return defs ->
        defs != null
            && defs.size() == expectedNames.size()
            && defs.stream().map(PipelineOutputDefinition::getName).toList().equals(expectedNames);
  }
}
