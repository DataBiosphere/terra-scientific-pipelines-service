package bio.terra.pipelines.common.utils;

import com.fasterxml.jackson.databind.ObjectMapper;

public enum PipelineInputTypesEnum {
  STRING("string") {
    @Override
    public String cast(Object value) {
      return objectMapper.convertValue(value, String.class);
    }
  },
  INTEGER("integer") {
    @Override
    public Integer cast(Object value) {
      // objectMapper will accept a float and convert it to an integer; we don't want this
      if (value instanceof Float || value instanceof Double) {
        throw new IllegalArgumentException("Integer input must be an integer");
      }
      return objectMapper.convertValue(value, Integer.class);
    }
  },
  ARRAY_STRING("array_string") {
    @Override
    public String[] cast(Object value) {
      return objectMapper.convertValue(value, String[].class);
    }
  };
  ObjectMapper objectMapper = new ObjectMapper();

  public abstract <T> T cast(Object value);

  private final String value;

  PipelineInputTypesEnum(String value) {
    this.value = value;
  }
}
