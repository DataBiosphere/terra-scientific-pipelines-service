package bio.terra.pipelines.dependencies.stairway;

public enum JobMapKeys {
  // parameters for all flight types
  DESCRIPTION("description"),
  USER_ID("user_id"),
  PIPELINE_ID("pipeline_id"),
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

  public static boolean isRequiredKey(String keyName) {
    return keyName.equals(JobMapKeys.DESCRIPTION.getKeyName())
        || keyName.equals(JobMapKeys.USER_ID.getKeyName())
        || keyName.equals(JobMapKeys.PIPELINE_ID.getKeyName());
  }
}
