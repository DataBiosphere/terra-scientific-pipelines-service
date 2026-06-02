package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import java.util.Objects;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Getter;
import lombok.ToString;

/**
 * Immutable domain model representing pipeline quota configuration. This is derived from the YAML
 * pipeline configuration and represents the semantic quota definition.
 */
@Getter
@ToString
@AllArgsConstructor
@Builder(toBuilder = true)
public class PipelineQuota {
  private final PipelinesEnum pipelineName;
  private final Integer defaultQuota;
  private final Integer minQuotaConsumed;
  private final QuotaUnitsEnum quotaUnits;

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    PipelineQuota that = (PipelineQuota) o;
    return Objects.equals(pipelineName, that.pipelineName)
        && Objects.equals(defaultQuota, that.defaultQuota)
        && Objects.equals(minQuotaConsumed, that.minQuotaConsumed)
        && quotaUnits == that.quotaUnits;
  }

  @Override
  public int hashCode() {
    return Objects.hash(pipelineName, defaultQuota, minQuotaConsumed, quotaUnits);
  }
}
