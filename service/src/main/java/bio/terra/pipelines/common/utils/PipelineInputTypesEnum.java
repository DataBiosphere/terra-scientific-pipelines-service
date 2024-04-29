package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.ValidationException;
import com.fasterxml.jackson.databind.ObjectMapper;

public enum PipelineInputTypesEnum {
  STRING() {
    @Override
    public String cast(String fieldName, Object value) {
      validateNotNull(fieldName, value);
      String stringValue = objectMapper.convertValue(value, String.class);
      if (stringValue.isBlank()) {
        throw new ValidationException(String.format("%s must not be empty", fieldName));
      }
      return stringValue;
    }
  },
  INTEGER() {
    @Override
    public Integer cast(String fieldName, Object value) {
      validateNotNull(fieldName, value);
      // objectMapper will accept a float and convert it to an integer; we don't want this
      if (value instanceof Float || value instanceof Double) {
        throw new IllegalArgumentException(
            "Integer input must be an integer"); // this error message gets overwritten
      }
      return objectMapper.convertValue(value, Integer.class);
    }
  },
  VCF() {
    @Override
    public String cast(String fieldName, Object value) {
      validateNotNull(fieldName, value);
      String stringValue = objectMapper.convertValue(value, String.class);
      if (!stringValue.endsWith(".vcf.gz")) {
        throw new ValidationException(
            String.format("%s must be a path to a VCF file ending in .vcf.gz", fieldName));
      }
      return stringValue;
    }
  },
  ARRAY_STRING() {
    @Override
    public String[] cast(String fieldName, Object value) {
      validateNotNull(fieldName, value);
      String[] stringArray = objectMapper.convertValue(value, String[].class);
      for (String stringValue : stringArray) {
        try {
          STRING.cast(fieldName, stringValue);
        } catch (ValidationException e) {
          throw new ValidationException(
              String.format("%s must be an array of non-null strings", fieldName));
        }
      }
      return stringArray;
    }
  },
  ARRAY_VCF() {
    @Override
    public String[] cast(String fieldName, Object value) {
      validateNotNull(fieldName, value);
      String[] stringArray = objectMapper.convertValue(value, String[].class);
      for (String stringValue : stringArray) {
        try {
          VCF.cast(fieldName, stringValue);
        } catch (ValidationException e) {
          throw new ValidationException(
              String.format(
                  "%s must be an array of non-null paths to VCF files ending in .vcf.gz",
                  fieldName));
        }
      }
      return stringArray;
    }
  };
  final ObjectMapper objectMapper = new ObjectMapper();

  public abstract <T> T cast(String fieldName, Object value);

  private static void validateNotNull(String fieldName, Object value) {
    if (value == null) {
      throw new ValidationException(String.format("%s must not be null", fieldName));
    }
  }
}
