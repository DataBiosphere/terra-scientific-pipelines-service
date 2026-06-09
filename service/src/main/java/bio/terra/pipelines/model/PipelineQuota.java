package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.ToString;

/**
 * Immutable domain model representing pipeline quota configuration. This is derived from the YAML
 * pipeline configuration and represents the semantic quota definition.
 */
@Getter
@ToString
@EqualsAndHashCode
@AllArgsConstructor
@Builder(toBuilder = true)
public class PipelineQuota {
  private final PipelinesEnum pipelineName;
  private final Integer defaultQuota;
  private final Integer minQuotaConsumed;
  private final QuotaUnitsEnum quotaUnits;
}
