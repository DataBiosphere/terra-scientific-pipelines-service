package bio.terra.pipelines.common.utils;

public enum PipelinesEnum {
  IMPUTATION("imputation");

  private final String value;

  PipelinesEnum(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
// public class PipelineIds {
//  public static final String IMPUTATION = "imputation";
//  public static final List<String> ALL_PIPELINES = List.of(IMPUTATION);
//
//  public static boolean pipelineExists(String pipelineId) {
//    return ALL_PIPELINES.contains(pipelineId);
//  }
//
//  private PipelineIds() {}
// }
