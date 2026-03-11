package bio.terra.pipelines.common.utils;

public enum DataDeliveryStatusEnum {
  RUNNING,
  SUCCEEDED,
  FAILED;

  public boolean isSuccess() {
    return this == SUCCEEDED;
  }

  public boolean isCompleted() {
    return this == SUCCEEDED || this == FAILED;
  }
}
