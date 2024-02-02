package bio.terra.pipelines.dependencies.stairway;

import java.util.Arrays;
import java.util.List;

public enum JobMapKeys {
  // parameters for all flight types
  DESCRIPTION("description"),
  USER_ID("user_id"),
  PIPELINE_NAME("pipeline_name"),
  STATUS_CODE("status_code"),
  RESPONSE("response"), // used for synchronous jobs
  RESULT_PATH("result_path"); // used for asynchronous jobs

  private final String keyName;

  JobMapKeys(String keyName) {
    this.keyName = keyName;
  }

  public String getKeyName() {
    return keyName;
  }

  public static List<String> getRequiredKeys() {
    return Arrays.asList(JobMapKeys.USER_ID.getKeyName(), JobMapKeys.PIPELINE_NAME.getKeyName());
  }
}
