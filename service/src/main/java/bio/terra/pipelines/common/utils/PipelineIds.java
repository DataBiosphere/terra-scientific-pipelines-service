package bio.terra.pipelines.common.utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PipelineIds {
  public static final String IMPUTATION = "imputation";
  public static final List<String> ALL_PIPELINES =
      Collections.unmodifiableList(Arrays.asList(IMPUTATION));

  public static boolean pipelineExists(String pipelineId) {
    return ALL_PIPELINES.contains(pipelineId);
  }

  private PipelineIds() {}
}
