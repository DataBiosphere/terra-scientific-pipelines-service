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
        throw new ValidationException(String.format("%s must be a string", fieldName));
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
          throw new ValidationException(String.format("%s must be an integer", fieldName));
        }
      } else {
        throw new ValidationException(String.format("%s must be an integer", fieldName));
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
              String.format("%s must be a path to a VCF file ending in .vcf.gz", fieldName));
        }
        return stringValue.trim();
      } else {
        throw new ValidationException(String.format("%s must be a string", fieldName));
      }
    }
  },
  STRING_ARRAY {
    @Override
    public List<String> cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof List listValue) {
        try {
          List stringList =
              listValue.stream().map(itemValue -> STRING.cast(fieldName, itemValue)).toList();
          validateNotEmptyList(fieldName, stringList);
          return stringList;
        } catch (ValidationException e) {
          throw new ValidationException(String.format("%s must be an array of strings", fieldName));
        }
      } else {
        throw new ValidationException(String.format("%s must be an array of strings", fieldName));
      }
    }
  },
  VCF_ARRAY {
    @Override
    public List<String> cast(String fieldName, Object value) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof List listValue) {
        try {
          List stringList =
              listValue.stream().map(itemValue -> VCF.cast(fieldName, itemValue)).toList();
          validateNotEmptyList(fieldName, stringList);
          return stringList;
        } catch (ValidationException e) {
          throw new ValidationException(
              String.format(
                  "%s must be an array of paths to VCF files ending in .vcf.gz", fieldName));
        }
      } else {
        throw new ValidationException(
            String.format(
                "%s must be an array of paths to VCF files ending in .vcf.gz", fieldName));
      }
    }
  };

  public abstract Object cast(String fieldName, Object value);

  private static void validateNotNullOrEmpty(String fieldName, Object value) {
    if (value == null) {
      throw new ValidationException(String.format("%s must not be null", fieldName));
    }
    if (value instanceof String stringValue && stringValue.isBlank()) {
      throw new ValidationException(String.format("%s must not be empty", fieldName));
    }
  }

  private static void validateNotEmptyList(String fieldName, List<?> listValue) {
    if (listValue.isEmpty()) {
      // note this error message will be overwritten to a more specific message in the catch block
      throw new ValidationException(String.format("%s must not be an empty list", fieldName));
    }
  }
}
