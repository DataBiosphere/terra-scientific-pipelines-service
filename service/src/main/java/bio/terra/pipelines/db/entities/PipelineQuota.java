package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import jakarta.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(
    name = "pipeline_quotas",
    uniqueConstraints = {@UniqueConstraint(columnNames = {"pipeline_name"})})
public class PipelineQuota {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_name", nullable = false)
  private PipelinesEnum pipelineName;

  @Column(name = "default_quota")
  private int defaultQuota;

  @Column(name = "min_quota_consumed")
  private int minQuotaConsumed;

  public PipelineQuota(PipelinesEnum pipelineName, int defaultQuota, int minQuotaConsumed) {
    this.pipelineName = pipelineName;
    this.defaultQuota = defaultQuota;
    this.minQuotaConsumed = minQuotaConsumed;
  }
}
