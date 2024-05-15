package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.ValidationException;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.List;
import java.util.Optional;

public enum PipelineInputTypesEnum {
  STRING {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      if (value instanceof String stringValue) {
        String trimmedString = stringValue.trim();
        if (!(trimmedString.isBlank())) {
          return (T) trimmedString;
        }
      }
      return null;
    }

    @Override
    public Optional<String> validate(String fieldName, Object value) {
      if (cast(fieldName, value, new TypeReference<String>() {}) == null) {
        return Optional.of("%s must be a string".formatted(fieldName));
      }
      return Optional.empty();
    }
  },
  INTEGER {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      if (value instanceof Integer integerValue) {
        return (T) integerValue;
      } else if (value instanceof String stringValue) {
        try {
          return (T) Integer.valueOf(stringValue);
        } catch (NumberFormatException e) {
          return null;
        }
      }
      return null;
    }

    @Override
    public Optional<String> validate(String fieldName, Object value) {
      if (cast(fieldName, value, new TypeReference<Integer>() {}) == null) {
        return Optional.of("%s must be an integer".formatted(fieldName));
      }
      return Optional.empty();
    }
  },
  VCF {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      return (T) STRING.cast(fieldName, value, new TypeReference<>() {});
    }

    @Override
    public Optional<String> validate(String fieldName, Object value) {
      String stringCastValue = STRING.cast(fieldName, value, new TypeReference<>() {});
      if (stringCastValue == null || !(stringCastValue.endsWith(".vcf.gz"))) {
        return Optional.of(
            "%s must be a path to a VCF file ending in .vcf.gz".formatted(fieldName));
      }
      return Optional.empty();
    }
  },
  STRING_ARRAY {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      List<?> listValue = processListValue(value);
      if (listValue == null) {
        return null;
      }
      List<String> stringList =
          listValue.stream()
              .map(itemValue -> STRING.cast(fieldName, itemValue, new TypeReference<String>() {}))
              .toList();
      // if any of the items in the cast list are null, it means the item was not cast-able as a
      // string, and we consider the entire list to be un-cast-able as a string array
      if (stringList.contains(null)) {
        return null;
      }
      return (T) stringList;
    }

    @Override
    public Optional<String> validate(String fieldName, Object value) {
      Optional<String> stringArrayErrorMessage =
          Optional.of("%s must be an array of strings".formatted(fieldName));
      if (PipelineInputTypesEnum.valueIsNullOrEmpty(value)) {
        return stringArrayErrorMessage;
      }

      List<String> listValue = cast(fieldName, value, new TypeReference<>() {});
      if (listValue == null || listValue.isEmpty()) {
        return stringArrayErrorMessage;
      }

      return Optional.empty();
    }
  },
  VCF_ARRAY {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      List<?> listValue = processListValue(value);
      if (listValue == null) {
        return null;
      }
      List<String> stringList =
          listValue.stream()
              .map(itemValue -> VCF.cast(fieldName, itemValue, new TypeReference<String>() {}))
              .toList();
      if (stringList.contains(null)) {
        return null;
      }

      return (T) stringList;
    }

    @Override
    public Optional<String> validate(String fieldName, Object value) {
      Optional<String> vcfArrayErrorMessage =
          Optional.of(
              "%s must be an array of paths to VCF files ending in .vcf.gz".formatted(fieldName));
      if (PipelineInputTypesEnum.valueIsNullOrEmpty(value)) {
        return vcfArrayErrorMessage;
      }

      List<?> listValue = processListValue(value);
      if (listValue == null || listValue.isEmpty()) {
        return vcfArrayErrorMessage;
      }

      // validate that all the items in the list are VCFs
      List<Optional<String>> validationMessages =
          listValue.stream().map(itemValue -> VCF.validate(fieldName, itemValue)).toList();
      if (!List.of(Optional.empty()).containsAll(validationMessages)) {
        return vcfArrayErrorMessage;
      }

      return Optional.empty();
    }
  };

  /**
   * Cast the value to the correct type. We assume the value has been validated before being cast.
   * If the value is cast-able, return the cast value. If the value is not cast-able, return null.
   *
   * @param fieldName - the name of the field being cast (used to construct error messages)
   * @param value - the value to cast
   * @param typeReference - a new instance of TypeReference
   */
  public abstract <T> T cast(String fieldName, Object value, TypeReference<T> typeReference);

  /**
   * Validate that the value is of the correct type. If the value passes validation, return null. If
   * the value fails validation, return a descriptive error message in a String.
   *
   * @param fieldName - the name of the field being validated (used to construct error messages)
   * @param value - the value to validate
   */
  public abstract Optional<String> validate(String fieldName, Object value)
      throws ValidationException;

  private static boolean valueIsNullOrEmpty(Object value) {
    if (value == null) {
      return true;
    }
    return value instanceof String stringValue && stringValue.isBlank();
  }

  private static List<String> convertStringToList(String stringyList) {
    ObjectMapper objectMapper = new ObjectMapper();
    try {
      return List.of(objectMapper.readValue(stringyList, String[].class));
    } catch (JsonProcessingException e) {
      return null;
    }
  }

  private static List<?> processListValue(Object value) {
    if (value instanceof String stringValue) {
      return convertStringToList(stringValue);
    } else if (value instanceof List<?> listValue) {
      return listValue;
    } else {
      return null;
    }
  }
}
