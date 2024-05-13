package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.ValidationException;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.List;

public enum PipelineInputTypesEnum {
  STRING {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof String stringValue) {
        return (T) stringValue.trim();
      } else {
        throw new ValidationException("%s must be a string".formatted(fieldName));
      }
    }
  },
  INTEGER {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof Integer integerValue) {
        return (T) integerValue;
      } else if (value instanceof String stringValue) {
        try {
          return (T) Integer.valueOf(stringValue);
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
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      validateNotNullOrEmpty(fieldName, value);
      String stringValue = "";
      try {
        stringValue = STRING.cast(fieldName, value, new TypeReference<>() {});
      } catch (ValidationException e) {
        throw new ValidationException("%s must be a string".formatted(fieldName));
      }
      if (!stringValue.endsWith(".vcf.gz")) {
        throw new ValidationException(
            "%s must be a path to a VCF file ending in .vcf.gz".formatted(fieldName));
      }
      return (T) stringValue;
    }
  },
  STRING_ARRAY {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof List<?> listValue) {
        try {
          List<?> stringList =
              listValue.stream()
                  .map(itemValue -> STRING.cast(fieldName, itemValue, new TypeReference<>() {}))
                  .toList();
          validateNotEmptyList(fieldName, stringList);
          return (T) stringList;
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
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      validateNotNullOrEmpty(fieldName, value);
      if (value instanceof List<?> listValue) {
        try {
          List<?> stringList =
              listValue.stream()
                  .map(itemValue -> VCF.cast(fieldName, itemValue, new TypeReference<>() {}))
                  .toList();
          validateNotEmptyList(fieldName, stringList);
          // Note: I wrapped this in List.copyOf to the unit test pass (with the extra check I
          // added)
          // Without this, the list it returns is a java.util.ImmutableCollections$ListN and the
          // test expects a
          // java.util.ImmutableCollections$List12.
          // I don't really know the difference and this might not really be necessary but just
          // leaving breadcrumbs.
          return (T) List.copyOf(stringList);
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

  public abstract <T> T cast(String fieldName, Object value, TypeReference<T> typeReference);

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
