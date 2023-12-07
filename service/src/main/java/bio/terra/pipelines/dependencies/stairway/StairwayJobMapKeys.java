package bio.terra.pipelines.dependencies.stairway;

public enum StairwayJobMapKeys {
  // parameters for all flight types
  DESCRIPTION("description"),
  REQUEST("request"),
  RESPONSE("response"),
  STATUS_CODE("status_code"),

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
        || keyName.equals(StairwayJobMapKeys.REQUEST.getKeyName());
  }
}
