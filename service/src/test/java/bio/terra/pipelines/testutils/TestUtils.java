package bio.terra.pipelines.testutils;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.*;
import bio.terra.rawls.model.MethodConfiguration;
import bio.terra.rawls.model.MethodRepoMethod;
import java.util.*;

/** A collection of utilities and constants useful for tests. */
public class TestUtils {
  // Pipelines test constants
  public static final PipelinesEnum TEST_PIPELINE_1_IMPUTATION_ENUM =
      PipelinesEnum.ARRAY_IMPUTATION;

  public static final Long TEST_PIPELINE_ID_1 = 1L;
  public static final int TEST_PIPELINE_VERSION_1 = 0;
  public static final String TEST_PIPELINE_DISPLAY_NAME_1 =
      "Test Pipeline Name"; // this matches the job pre-populated in the db for tests
  public static final String TEST_PIPELINE_DESCRIPTION_1 = "Test Pipeline Description";
  public static final String TEST_PIPELINE_TYPE_1 = "imputation1";
  public static final String TEST_WDL_URL_1 = "http://nowhere1";
  public static final String TEST_TOOL_NAME_1 = "methodName1";
  public static final String TEST_TOOL_VERSION_1 = "0.1.12";

  public static final int TEST_PIPELINE_VERSION_2 = 1;
  public static final String TEST_PIPELINE_DISPLAY_NAME_2 = "Test Pipeline Name Two";
  public static final String TEST_PIPELINE_DESCRIPTION_2 = "Test Pipeline Description Two";
  public static final String TEST_PIPELINE_TYPE_2 = "imputation2";
  public static final String TEST_WDL_URL_2 = "http://nowhere2";
  public static final String TEST_TOOL_NAME_2 = "methodName2";
  public static final String TEST_TOOL_VERSION_2 = "1.1.12";
  public static final UUID CONTROL_WORKSPACE_ID =
      UUID.fromString("fafafafa-fafa-fafa-fafa-fafafafafafa");
  public static final String CONTROL_WORKSPACE_BILLING_PROJECT = "testTerraProject";
  public static final String CONTROL_WORKSPACE_NAME = "testTerraWorkspaceName";
  public static final String CONTROL_WORKSPACE_CONTAINER_NAME =
      "fc-secure-%s".formatted(CONTROL_WORKSPACE_ID);
  public static final String GCP_STORAGE_PROTOCOL = "gs://";
  public static final String CONTROL_WORKSPACE_GOOGLE_PROJECT = "testGoogleProject";
  public static final Map<String, String> TEST_PIPELINE_OUTPUTS =
      new HashMap(
          Map.of(
              "testFileOutputKey",
              "gs://fc-secure-%s/testFileOutputValue".formatted(CONTROL_WORKSPACE_ID)));

  public static final List<PipelineInputDefinition> TEST_PIPELINE_INPUTS_DEFINITION_LIST =
      new ArrayList<>(
          List.of(
              new PipelineInputDefinition(
                  3L,
                  "testRequiredStringInput",
                  "test_required_string_input",
                  PipelineVariableTypesEnum.STRING,
                  null,
                  true,
                  true,
                  false,
                  null),
              new PipelineInputDefinition(
                  3L,
                  "testOptionalStringInput",
                  "test_optional_string_input",
                  PipelineVariableTypesEnum.STRING,
                  null,
                  false,
                  true,
                  false,
                  "testDefaultValue"),
              new PipelineInputDefinition(
                  3L,
                  "testRequiredIntInput",
                  "test_required_int_input",
                  PipelineVariableTypesEnum.INTEGER,
                  null,
                  true,
                  true,
                  false,
                  null),
              new PipelineInputDefinition(
                  3L,
                  "testOptionalIntInput",
                  "test_optional_int_input",
                  PipelineVariableTypesEnum.INTEGER,
                  null,
                  false,
                  true,
                  false,
                  "42"),
              new PipelineInputDefinition(
                  3L,
                  "testServiceProvidedInput",
                  "test_service_provided_input",
                  PipelineVariableTypesEnum.STRING,
                  null,
                  true,
                  false,
                  false,
                  "testServiceProvidedDefaultValue"),
              new PipelineInputDefinition(
                  3L,
                  "testRequiredVcfInput",
                  "test_required_vcf_input",
                  PipelineVariableTypesEnum.FILE,
                  ".vcf.gz",
                  true,
                  true,
                  false,
                  null)));

  public static final List<PipelineOutputDefinition> TEST_PIPELINE_OUTPUTS_DEFINITION_LIST =
      new ArrayList<>(
          List.of(
              new PipelineOutputDefinition(
                  3L, "outputName", "output_name", PipelineVariableTypesEnum.FILE)));

  public static final Pipeline TEST_PIPELINE_1 =
      new Pipeline(
          PipelinesEnum.ARRAY_IMPUTATION,
          TEST_PIPELINE_VERSION_1,
          TEST_PIPELINE_DISPLAY_NAME_1,
          TEST_PIPELINE_DESCRIPTION_1,
          TEST_PIPELINE_TYPE_1,
          TEST_WDL_URL_1,
          TEST_TOOL_NAME_1,
          TEST_TOOL_VERSION_1,
          CONTROL_WORKSPACE_BILLING_PROJECT,
          CONTROL_WORKSPACE_NAME,
          CONTROL_WORKSPACE_CONTAINER_NAME,
          CONTROL_WORKSPACE_GOOGLE_PROJECT,
          TEST_PIPELINE_INPUTS_DEFINITION_LIST,
          TEST_PIPELINE_OUTPUTS_DEFINITION_LIST);
  public static final Pipeline TEST_PIPELINE_2 =
      new Pipeline(
          PipelinesEnum.ARRAY_IMPUTATION,
          TEST_PIPELINE_VERSION_2,
          TEST_PIPELINE_DISPLAY_NAME_2,
          TEST_PIPELINE_DESCRIPTION_2,
          TEST_PIPELINE_TYPE_2,
          TEST_WDL_URL_2,
          TEST_TOOL_NAME_2,
          TEST_TOOL_VERSION_2,
          CONTROL_WORKSPACE_BILLING_PROJECT,
          CONTROL_WORKSPACE_NAME,
          CONTROL_WORKSPACE_CONTAINER_NAME,
          CONTROL_WORKSPACE_GOOGLE_PROJECT,
          TEST_PIPELINE_INPUTS_DEFINITION_LIST,
          TEST_PIPELINE_OUTPUTS_DEFINITION_LIST);

  public static Pipeline createTestPipelineWithId() {
    Pipeline pipeline =
        new Pipeline(
            TestUtils.TEST_PIPELINE_1.getName(),
            TestUtils.TEST_PIPELINE_1.getVersion(),
            TestUtils.TEST_PIPELINE_1.getDisplayName(),
            TestUtils.TEST_PIPELINE_1.getDescription(),
            TestUtils.TEST_PIPELINE_1.getPipelineType(),
            TestUtils.TEST_PIPELINE_1.getWdlUrl(),
            TestUtils.TEST_PIPELINE_1.getToolName(),
            TestUtils.TEST_PIPELINE_1.getToolVersion(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceBillingProject(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceName(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceStorageContainerName(),
            TestUtils.TEST_PIPELINE_1.getWorkspaceGoogleProject(),
            TestUtils.TEST_PIPELINE_1.getPipelineInputDefinitions(),
            TestUtils.TEST_PIPELINE_1.getPipelineOutputDefinitions());
    pipeline.setId(1L);
    return pipeline;
  }

  public static final PipelineQuota TEST_PIPELINE_QUOTA_1 =
      new PipelineQuota(PipelinesEnum.ARRAY_IMPUTATION, 100, 10, QuotaUnitsEnum.SAMPLES);

  public static final String TEST_USER_ID_1 =
      "testUser"; // this matches the job pre-populated in the db for tests
  public static final String TEST_USER_ID_2 = "testUser2";

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
    return createNewPipelineRunWithJobIdAndUser(jobId, TEST_USER_ID_1);
  }

  public static PipelineRun createNewPipelineRunWithJobIdAndUser(UUID jobId, String userId) {
    return new PipelineRun(
        jobId,
        userId,
        TEST_PIPELINE_ID_1,
        TEST_TOOL_VERSION_1,
        CONTROL_WORKSPACE_BILLING_PROJECT,
        CONTROL_WORKSPACE_NAME,
        CONTROL_WORKSPACE_CONTAINER_NAME,
        CONTROL_WORKSPACE_GOOGLE_PROJECT,
        CommonPipelineRunStatusEnum.PREPARING,
        "test description");
  }

  public static String buildTestResultUrl(String jobId) {
    return "https://%s/api/pipelineruns/v1/result/%s".formatted(TEST_DOMAIN, jobId);
  }
}
