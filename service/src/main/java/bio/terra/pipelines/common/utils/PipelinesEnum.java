package bio.terra.pipelines.common.utils;

public enum PipelinesEnum {
  IMPUTATION_MINIMAC4("imputation_minimac4");

  private final String value;

  PipelinesEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
