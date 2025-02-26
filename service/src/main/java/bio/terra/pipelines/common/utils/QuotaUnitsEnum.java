package bio.terra.pipelines.common.utils;

public enum QuotaUnitsEnum {
  SAMPLES("samples");

  private final String value;

  QuotaUnitsEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
