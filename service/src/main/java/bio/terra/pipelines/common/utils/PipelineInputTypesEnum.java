package bio.terra.pipelines.common.utils;

public enum PipelineInputTypesEnum {
  STRING("String") {
    @Override
    public String cast(Object value) {
      return String.valueOf(value);
    }
  },
  INTEGER("Integer") {
    @Override
    public Integer cast(Object value) {
      return Integer.parseInt(String.valueOf(value));
    }
  };

  public abstract <T> T cast(Object value);

  private final String value;

  PipelineInputTypesEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
