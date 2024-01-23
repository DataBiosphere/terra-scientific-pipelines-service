package bio.terra.pipelines.testutils;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.UUID;

/** A collection of utilities and constants useful for tests. */
public class TestUtils {
  // Pipelines test constants
  public static final PipelinesEnum TEST_PIPELINE_1_ENUM = PipelinesEnum.IMPUTATION;
  public static final String TEST_PIPELINE_ID_1 =
      TEST_PIPELINE_1_ENUM
          .getValue(); // this matches the job pre-populated in the db for tests in that it is in
  // the imputation_jobs table
  public static final String TEST_PIPELINE_VERSION_1 =
      "testVersion"; // this matches the job pre-populated in the db for tests
  public static final String TEST_PIPELINE_NAME_1 =
      "Test Pipeline Name"; // this matches the job pre-populated in the db for tests
  public static final String TEST_PIPELINE_DESCRIPTION_1 = "Test Pipeline Description";
  public static final String TEST_PIPELINE_ID_2 = "testPipeline2";
  public static final String TEST_PIPELINE_VERSION_2 = "testVersion2";
  public static final String TEST_PIPELINE_NAME_2 = "Test Pipeline Name Two";
  public static final String TEST_PIPELINE_DESCRIPTION_2 = "Test Pipeline Description Two";
  public static final Pipeline TEST_PIPELINE_1 =
      new Pipeline(
          TEST_PIPELINE_ID_1,
          TEST_PIPELINE_VERSION_1,
          TEST_PIPELINE_NAME_1,
          TEST_PIPELINE_DESCRIPTION_1);
  public static final Pipeline TEST_PIPELINE_2 =
      new Pipeline(
          TEST_PIPELINE_ID_2,
          TEST_PIPELINE_VERSION_2,
          TEST_PIPELINE_NAME_2,
          TEST_PIPELINE_DESCRIPTION_2);

  public static final String TEST_USER_ID_1 =
      "testUser"; // this matches the job pre-populated in the db for tests
  public static final String TEST_USER_ID_2 = "testUser2";

  public static final UUID TEST_EXISTING_UUID =
      // this matches the job pre-populated in the db for tests
      UUID.fromString("deadbeef-dead-beef-deaf-beefdeadbeef");

  public static final UUID TEST_NEW_UUID = UUID.fromString("deadbeef-dead-beef-aaaa-beefdeadbeef");

  public static final String TEST_STATUS = "TEST STATUS";

  public static final Object TEST_PIPELINE_INPUTS =
      new LinkedHashMap<>(Map.of("first_key", "first_value"));
}
