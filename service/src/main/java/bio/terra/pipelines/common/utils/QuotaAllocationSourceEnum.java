package bio.terra.pipelines.common.utils;

public enum QuotaAllocationSourceEnum {
  DEFAULT_FREE("Default_Free"),
  FNF("FNF"),
  DEVELOPMENT_FREE("Development_Free"),
  PRODUCTION_FREE("Production_Free"),
  PRODUCTION_CLIENT_TESTING("Production_Client_Testing"),
  PAID_EXTERNAL("Paid_External"),
  PAID_INTERNAL("Paid_Internal");

  private final String value;

  QuotaAllocationSourceEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
