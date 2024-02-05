package bio.terra.pipelines.app.common;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import io.micrometer.core.instrument.Metrics;
import lombok.AccessLevel;
import lombok.NoArgsConstructor;

/**
 * Class used to interact with micrometer metrics that are exported to /actuator/prometheus. This
 * was mostly taken from <a
 * href="https://github.com/DataBiosphere/terra-billing-profile-manager/blob/main/service/src/main/java/bio/terra/profile/app/common/MetricUtils.java">bpm
 * metrics</a>
 */
@NoArgsConstructor(access = AccessLevel.PRIVATE)
public class MetricsUtils {
  private static final String NAMESPACE = "tsps";
  private static final String PIPELINE_TAG = "pipeline";

  /**
   * increments metrics counter for a tsps pipeline that has been run
   *
   * @param pipelineName - name of pipeline that was run e.g. "imputation"
   */
  public static void incrementPipelineRun(PipelinesEnum pipelineName) {
    Metrics.globalRegistry
        .counter(
            String.format("%s.pipeline.run.count", NAMESPACE),
            PIPELINE_TAG,
            pipelineName.getValue())
        .increment();
  }
}
