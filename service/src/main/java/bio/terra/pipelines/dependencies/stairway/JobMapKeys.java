package bio.terra.pipelines.dependencies.stairway;

public enum JobMapKeys {
  // parameters for all flight types
  DESCRIPTION("description"),
  REQUEST("request"),
  USER_ID("user_id"),
  PIPELINE_ID("pipeline_id"),
  RESPONSE("response"), // TODO what will this actually be used for?
  STATUS_CODE("status_code"),
  RESULT_PATH("result_path"),

  // parameter for the job
  FLIGHT_CLASS("flight_class");

  private final String keyName;

  JobMapKeys(String keyName) {
    this.keyName = keyName;
  }

  public String getKeyName() {
    return keyName;
  }

  // TODO use this
  public static boolean isRequiredKey(String keyName) {
    return keyName.equals(JobMapKeys.DESCRIPTION.getKeyName())
        || keyName.equals(JobMapKeys.REQUEST.getKeyName())
        || keyName.equals(JobMapKeys.USER_ID.getKeyName())
        || keyName.equals(JobMapKeys.PIPELINE_ID.getKeyName());
  }
}
