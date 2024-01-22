package bio.terra.pipelines.dependencies.stairway;

public enum StairwayJobMapKeys {
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

  StairwayJobMapKeys(String keyName) {
    this.keyName = keyName;
  }

  public String getKeyName() {
    return keyName;
  }

  public static boolean isRequiredKey(String keyName) {
    return keyName.equals(StairwayJobMapKeys.DESCRIPTION.getKeyName())
        || keyName.equals(StairwayJobMapKeys.REQUEST.getKeyName())
        || keyName.equals(StairwayJobMapKeys.USER_ID.getKeyName())
        || keyName.equals(StairwayJobMapKeys.PIPELINE_ID.getKeyName());
  }
}
