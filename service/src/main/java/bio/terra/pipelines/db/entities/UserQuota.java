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
    name = "user_quotas",
    uniqueConstraints = {@UniqueConstraint(columnNames = {"pipeline_name", "user_id"})})
public class UserQuota {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_name", nullable = false)
  private PipelinesEnum pipelineName;

  @Column(name = "user_id", nullable = false)
  private String userId;

  @Column(name = "quota")
  private int quota;

  @Column(name = "quota_consumed")
  private int quotaConsumed;

  public UserQuota(PipelinesEnum pipelineName, String userId, int quota, int quotaConsumed) {
    this.pipelineName = pipelineName;
    this.userId = userId;
    this.quota = quota;
    this.quotaConsumed = quotaConsumed;
  }
}
