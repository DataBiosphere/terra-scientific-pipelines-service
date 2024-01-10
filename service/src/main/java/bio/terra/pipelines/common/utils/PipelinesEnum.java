package bio.terra.pipelines.common.utils;

public enum PipelinesEnum {
  IMPUTATION("imputation");

  private final String value;

  PipelinesEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
