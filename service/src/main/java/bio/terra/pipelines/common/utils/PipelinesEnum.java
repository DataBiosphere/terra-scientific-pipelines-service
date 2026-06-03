package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.NotFoundException;
import lombok.Getter;

public enum PipelinesEnum {
  ARRAY_IMPUTATION("array_imputation", "arrayImputation"),
  LOW_PASS_IMPUTATION("low_pass_imputation", "lowPassImputation");

  @Getter private final String lowerCaseValue;
  @Getter private final String configKeyValue;

  PipelinesEnum(String lowerCaseValue, String configKeyValue) {
    this.lowerCaseValue = lowerCaseValue;
    this.configKeyValue = configKeyValue;
  }

  public static PipelinesEnum enumFromLowerCaseValue(String stringValue) {
    return PipelinesEnum.valueOf(stringValue.toUpperCase());
  }

  public static PipelinesEnum enumFromConfigKeyValue(String configKeyValue) {
    return switch (configKeyValue) {
      case "arrayImputation" -> PipelinesEnum.ARRAY_IMPUTATION;
      case "lowPassImputation" -> PipelinesEnum.LOW_PASS_IMPUTATION;
      default ->
          throw new NotFoundException(
              "Pipeline not found for configKeyValue %s".formatted(configKeyValue));
    };
  }
}
