package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.ValidationException;
import java.util.List;

public enum PipelineInputTypesEnum {
  STRING {
    @Override
    public String cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof String stringValue) {
        return stringValue.trim();
      } else {
        throw new ValidationException("%s must be a string".formatted(fieldName));
      }
    }
  },
  INTEGER {
    @Override
    public Integer cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof Integer integerValue) {
        return integerValue;
      } else if (value instanceof String stringValue) {
        try {
          return Integer.parseInt(stringValue);
        } catch (NumberFormatException e) {
          throw new ValidationException("%s must be an integer".formatted(fieldName));
        }
      } else {
        throw new ValidationException("%s must be an integer".formatted(fieldName));
      }
    }
  },
  VCF {
    @Override
    public String cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof String stringValue) {
        if (!stringValue.trim().endsWith(".vcf.gz")) {
          throw new ValidationException(
              "%s must be a path to a VCF file ending in .vcf.gz".formatted(fieldName));
        }
        return stringValue.trim();
      } else {
        throw new ValidationException("%s must be a string".formatted(fieldName));
      }
    }
  },
  STRING_ARRAY {
    @Override
    public List<String> cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof List<?> listValue) {
        try {
          List<String> stringList =
              listValue.stream()
                  .map(itemValue -> STRING.cast(fieldName, itemValue).toString())
                  .toList();
          validateNotEmptyList(fieldName, stringList);
          return stringList;
        } catch (ValidationException e) {
          throw new ValidationException("%s must be an array of strings".formatted(fieldName));
        }
      } else {
        throw new ValidationException("%s must be an array of strings".formatted(fieldName));
      }
    }
  },
  VCF_ARRAY {
    @Override
    public List<String> cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof List<?> listValue) {
        try {
          List<String> stringList =
              listValue.stream()
                  .map(itemValue -> VCF.cast(fieldName, itemValue).toString())
                  .toList();
          validateNotEmptyList(fieldName, stringList);
          return stringList;
        } catch (ValidationException e) {
          throw new ValidationException(
              "%s must be an array of paths to VCF files ending in .vcf.gz".formatted(fieldName));
        }
      } else {
        throw new ValidationException(
            "%s must be an array of paths to VCF files ending in .vcf.gz".formatted(fieldName));
      }
    }
  };

  public abstract Object cast(String fieldName, Object value);

  private static void validateNotNullOrEmpty(String fieldName, Object value) {
    if (value == null) {
      throw new ValidationException("%s must not be null".formatted(fieldName));
    }
    if (value instanceof String stringValue && stringValue.isBlank()) {
      throw new ValidationException("%s must not be empty".formatted(fieldName));
    }
  }

  private static void validateNotEmptyList(String fieldName, List<?> listValue) {
    if (listValue.isEmpty()) {
      // note this error message will be overwritten to a more specific message in the catch block
      throw new ValidationException("%s must not be an empty list".formatted(fieldName));
    }
  }
}
