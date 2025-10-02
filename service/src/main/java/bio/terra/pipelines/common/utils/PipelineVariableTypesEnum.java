package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.List;
import java.util.Objects;
import java.util.regex.Pattern;

public enum PipelineVariableTypesEnum {
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
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      String castValue = cast(fieldName, value, new TypeReference<String>() {});
      if (castValue == null) {
        return "%s must be a string".formatted(fieldName);
      }
      if (!VALID_STRING_PATTERN.matcher(castValue).matches()) {
        return "%s must only contain alphanumeric characters or the following symbols: -_.=\\/"
            .formatted(fieldName);
      }
      return null;
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
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      Integer castValue = cast(fieldName, value, new TypeReference<>() {});
      if (castValue == null) {
        return "%s must be an integer".formatted(fieldName);
      }
      return PipelineVariableTypesEnum.checkRange(
          Double.valueOf(castValue), pipelineInputDefinition);
    }
  },
  FLOAT {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      if (value instanceof Double doubleValue) {
        return (T) doubleValue;
      } else if (value instanceof Integer integerValue) {
        return (T) Double.valueOf(integerValue);
      } else if (value instanceof String stringValue) {
        try {
          return (T) Double.valueOf(stringValue);
        } catch (NumberFormatException e) {
          return null;
        }
      }
      return null;
    }

    @Override
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      Double castValue = cast(fieldName, value, new TypeReference<>() {});
      if (cast(fieldName, value, new TypeReference<Double>() {}) == null) {
        return "%s must be a float".formatted(fieldName);
      }
      return PipelineVariableTypesEnum.checkRange(castValue, pipelineInputDefinition);
    }
  },
  BOOLEAN {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      if (value instanceof Boolean booleanValue) {
        return (T) booleanValue;
      } else if (value instanceof String stringValue) {
        String trimmedString = stringValue.trim().toLowerCase();
        if (trimmedString.equals("true")) {
          return (T) Boolean.TRUE;
        } else if (trimmedString.equals("false")) {
          return (T) Boolean.FALSE;
        }
      }
      return null;
    }

    @Override
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      if (cast(fieldName, value, new TypeReference<Boolean>() {}) == null) {
        return "%s must be a boolean".formatted(fieldName);
      }
      return null;
    }
  },
  FILE {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      return (T) STRING.cast(fieldName, value, new TypeReference<>() {});
    }

    @Override
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      String fileSuffix = pipelineInputDefinition.getFileSuffix();
      String stringCastValue = STRING.cast(fieldName, value, new TypeReference<>() {});
      if (stringCastValue == null || !(stringCastValue.endsWith(fileSuffix))) {
        return "%s must be a path to a file ending in %s".formatted(fieldName, fileSuffix);
      }
      if (!VALID_STRING_PATTERN.matcher(stringCastValue).matches()) {
        return "%s must only contain alphanumeric characters or the following symbols: -_.=\\/"
            .formatted(fieldName);
      }
      return null;
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
      // if any item in the cast list is null, it means the item was not cast-able as a
      // string, and we consider the entire list to be un-cast-able as a string array
      if (stringList.contains(null)) {
        return null;
      }
      return (T) stringList;
    }

    @Override
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      String stringArrayErrorMessage = "%s must be an array of strings".formatted(fieldName);
      if (value == null) {
        return NOT_NULL_OR_EMPTY_ERROR_MESSAGE.formatted(fieldName);
      }

      List<String> listValue = cast(fieldName, value, new TypeReference<>() {});
      if (listValue == null) {
        return stringArrayErrorMessage;
      } else if (listValue.isEmpty()) {
        return NOT_NULL_OR_EMPTY_ERROR_MESSAGE.formatted(fieldName);
      }

      // validate each string in the array
      for (String item : listValue) {
        if (!VALID_STRING_PATTERN.matcher(item).matches()) {
          return "%s must only contain strings with alphanumeric characters or the following symbols: -_.=\\/"
              .formatted(fieldName);
        }
      }

      // no issues found
      return null;
    }
  },
  FILE_ARRAY {
    @Override
    public <T> T cast(String fieldName, Object value, TypeReference<T> typeReference) {
      List<?> listValue = processListValue(value);
      if (listValue == null) {
        return null;
      }
      List<String> stringList =
          listValue.stream()
              .map(itemValue -> FILE.cast(fieldName, itemValue, new TypeReference<String>() {}))
              .toList();
      // if any item in the cast list is null, it means the item was not cast-able as a
      // string, and we consider the entire list to be un-cast-able as a FILE array.
      // Note that the FILE cast does not include a suffix check.
      if (stringList.contains(null)) {
        return null;
      }

      return (T) stringList;
    }

    @Override
    public String validate(PipelineInputDefinition pipelineInputDefinition, Object value) {
      String fieldName = pipelineInputDefinition.getName();
      String fileSuffix = pipelineInputDefinition.getFileSuffix();
      String fileArrayErrorMessage =
          ("%s must be an array of paths to files ending in %s and containing only alphanumeric characters or the following symbols: -_.=\\/")
              .formatted(fieldName, fileSuffix);
      if (value == null) {
        return NOT_NULL_OR_EMPTY_ERROR_MESSAGE.formatted(fieldName);
      }

      List<?> listValue = processListValue(value);
      if (listValue == null) {
        return fileArrayErrorMessage;
      } else if (listValue.isEmpty()) {
        return NOT_NULL_OR_EMPTY_ERROR_MESSAGE.formatted(fieldName);
      }

      // validate that all the items in the list are FILEs with the correct suffix and only contain
      // valid characters
      List<String> validationMessages =
          listValue.stream()
              .map(itemValue -> FILE.validate(pipelineInputDefinition, itemValue))
              .toList();
      if (validationMessages.stream().anyMatch(Objects::nonNull)) {
        return fileArrayErrorMessage;
      }

      // no issues found
      return null;
    }
  };

  /**
   * Cast the value to the correct type. We assume the value has been validated before being cast.
   * If the value is cast-able, return the cast value. If the value is not cast-able, return null.
   *
   * @param fieldName - the name of the field being cast (used to construct error messages)
   * @param value - the value to cast
   * @param typeReference - TypeReference indicating the type to be returned
   */
  public abstract <T> T cast(String fieldName, Object value, TypeReference<T> typeReference);

  /**
   * Validate that the value is of the correct type. If the value passes validation, return null. If
   * the value fails validation, return a descriptive error message in a String.
   *
   * @param pipelineInputDefinition - the input definition that contains parameters for validation
   *     (fieldName, fileSuffix)
   * @param value - the value to validate
   */
  public abstract String validate(PipelineInputDefinition pipelineInputDefinition, Object value);

  private static final String NOT_NULL_OR_EMPTY_ERROR_MESSAGE = "%s must not be null or empty";
  // this regex only allows alphanumeric characters, dashes, underscores, periods, equal signs, and
  // forward and backward slashes
  private static final Pattern VALID_STRING_PATTERN = Pattern.compile("^[a-zA-Z0-9_.=\\\\/-]+$");

  @SuppressWarnings(
      "java:S1168") // Disable "Empty arrays and collections should be returned instead of null"
  private static List<String> convertStringToList(String stringyList) {
    ObjectMapper objectMapper = new ObjectMapper();
    try {
      return List.of(objectMapper.readValue(stringyList, String[].class));
    } catch (JsonProcessingException e) {
      return null;
    }
  }

  @SuppressWarnings(
      "java:S1168") // Disable "Empty arrays and collections should be returned instead of null"
  private static List<?> processListValue(Object value) {
    if (value instanceof String stringValue) {
      return convertStringToList(stringValue);
    } else if (value instanceof List<?> listValue) {
      return listValue;
    } else {
      return null;
    }
  }

  // method that checks if cast value is in between min and max, inclusive on both ends, of pipeline
  // input definition if either is set and returns an appropriate error message if not in between,
  // otherwise returns null
  private static String checkRange(
      Double castValue, PipelineInputDefinition pipelineInputDefinition) {
    Double min = pipelineInputDefinition.getMinValue();
    Double max = pipelineInputDefinition.getMaxValue();
    String fieldName = pipelineInputDefinition.getName();

    boolean belowMin = min != null && castValue.compareTo(min) < 0;
    boolean aboveMax = max != null && castValue.compareTo(max) > 0;

    if (!belowMin && !aboveMax) return null;

    if (min != null && max != null) {
      return "%s must be between %s and %s".formatted(fieldName, min, max);
    } else if (min != null) {
      return "%s must be at least %s".formatted(fieldName, min);
    } else {
      return "%s must be at most %s".formatted(fieldName, max);
    }
  }
}
