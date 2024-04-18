package bio.terra.pipelines.common.utils;

public enum PipelinesEnum {
  IMPUTATION_BEAGLE("imputation_beagle");

  private final String value;

  PipelinesEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
