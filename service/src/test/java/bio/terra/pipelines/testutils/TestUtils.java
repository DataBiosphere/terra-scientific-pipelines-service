package bio.terra.pipelines.testutils;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import java.util.*;

/** A collection of utilities and constants useful for tests. */
public class TestUtils {
  // Pipelines test constants
  public static final PipelinesEnum TEST_PIPELINE_1_IMPUTATION_ENUM =
      PipelinesEnum.IMPUTATION_BEAGLE;
  public static final String TEST_PIPELINE_NAME_1_IMPUTATION =
      TEST_PIPELINE_1_IMPUTATION_ENUM
          .getValue(); // this matches the job pre-populated in the db for tests in that it is in
  // the jobs table

  public static final Long TEST_PIPELINE_ID_1 = 3L;
  public static final String TEST_PIPELINE_VERSION_1 =
      "testVersion"; // this matches the job pre-populated in the db for tests
  public static final String TEST_PIPELINE_DISPLAY_NAME_1 =
      "Test Pipeline Name"; // this matches the job pre-populated in the db for tests
  public static final String TEST_PIPELINE_DESCRIPTION_1 = "Test Pipeline Description";
  public static final String TEST_PIPELINE_TYPE_1 = "imputation1";
  public static final String TEST_WDL_URL_1 = "http://nowhere1";
  public static final String TEST_WDL_METHOD_NAME_1 = "methodName1";
  public static final UUID TEST_WORKSPACE_ID_1 = UUID.randomUUID();
  public static final String TEST_PIPELINE_ID_2 = "testPipeline2";
  public static final String TEST_PIPELINE_VERSION_2 = "testVersion2";
  public static final String TEST_PIPELINE_DISPLAY_NAME_2 = "Test Pipeline Name Two";
  public static final String TEST_PIPELINE_DESCRIPTION_2 = "Test Pipeline Description Two";

  public static final String TEST_PIPELINE_TYPE_2 = "imputation2";
  public static final String TEST_WDL_URL_2 = "http://nowhere2";
  public static final String TEST_WDL_METHOD_NAME_2 = "methodName2";
  public static final UUID TEST_WORKSPACE_ID_2 = UUID.randomUUID();
  public static final UUID CONTROL_WORKSPACE_ID =
      UUID.fromString("fafafafa-fafa-fafa-fafa-fafafafafafa");
  public static final Map<String, String> TEST_PIPELINE_OUTPUTS =
      new HashMap(
          Map.of(
              "testFileOutputKey",
              "https://lz123.stuff/sc-%s/testFileOutputValue".formatted(CONTROL_WORKSPACE_ID)));

  public static final List<PipelineInputDefinition> TEST_PIPELINE_INPUTS_DEFINITION_LIST =
      new ArrayList<>(
          List.of(
              new PipelineInputDefinition(
                  3L,
                  "testRequiredStringInput",
                  "test_required_string_input",
                  PipelineInputTypesEnum.STRING,
                  true,
                  true,
                  null),
              new PipelineInputDefinition(
                  3L,
                  "testOptionalStringInput",
                  "test_optional_string_input",
                  PipelineInputTypesEnum.STRING,
                  false,
                  true,
                  "testDefaultValue"),
              new PipelineInputDefinition(
                  3L,
                  "testRequiredIntInput",
                  "test_required_int_input",
                  PipelineInputTypesEnum.INTEGER,
                  true,
                  true,
                  null),
              new PipelineInputDefinition(
                  3L,
                  "testOptionalIntInput",
                  "test_optional_int_input",
                  PipelineInputTypesEnum.INTEGER,
                  false,
                  true,
                  "42"),
              new PipelineInputDefinition(
                  3L,
                  "testServiceProvidedInput",
                  "test_service_provided_input",
                  PipelineInputTypesEnum.STRING,
                  true,
                  false,
                  "testServiceProvidedDefaultValue"),
              new PipelineInputDefinition(
                  3L,
                  "testRequiredVcfInput",
                  "test_required_vcf_input",
                  PipelineInputTypesEnum.VCF,
                  true,
                  true,
                  null)));
  ;
  public static final Pipeline TEST_PIPELINE_1 =
      new Pipeline(
          TEST_PIPELINE_NAME_1_IMPUTATION,
          TEST_PIPELINE_VERSION_1,
          TEST_PIPELINE_DISPLAY_NAME_1,
          TEST_PIPELINE_DESCRIPTION_1,
          TEST_PIPELINE_TYPE_1,
          TEST_WDL_URL_1,
          TEST_WDL_METHOD_NAME_1,
          TEST_WORKSPACE_ID_1,
          TEST_PIPELINE_INPUTS_DEFINITION_LIST);
  public static final Pipeline TEST_PIPELINE_2 =
      new Pipeline(
          TEST_PIPELINE_ID_2,
          TEST_PIPELINE_VERSION_2,
          TEST_PIPELINE_DISPLAY_NAME_2,
          TEST_PIPELINE_DESCRIPTION_2,
          TEST_PIPELINE_TYPE_2,
          TEST_WDL_URL_2,
          TEST_WDL_METHOD_NAME_2,
          TEST_WORKSPACE_ID_2,
          TEST_PIPELINE_INPUTS_DEFINITION_LIST);

  public static final String TEST_USER_ID_1 =
      "testUser"; // this matches the job pre-populated in the db for tests
  public static final String TEST_USER_ID_2 = "testUser2";

  public static final UUID TEST_NEW_UUID = UUID.fromString("deadbeef-dead-beef-aaaa-beefdeadbeef");

  public static final UUID TEST_NEW_UUID_2 =
      UUID.fromString("deadbeef-dead-beef-bbbb-beefdeadbeef");

  public static final Map<String, Object> TEST_PIPELINE_INPUTS =
      new HashMap<>(Map.of("first_key", "first_value"));

  public static final Map<String, Object> TEST_PIPELINE_INPUTS_IMPUTATION_BEAGLE =
      new HashMap<>(
          Map.of("multiSampleVcf", "fake/file.vcf.gz", "outputBasename", "fake_basename"));

  public static final String TEST_RESULT_URL = "https://some-teaspoons-domain.com/test/result/path";
  public static final String TEST_DOMAIN = "some-teaspoons-domain.com";
}
