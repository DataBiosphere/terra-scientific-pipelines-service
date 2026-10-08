package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import jakarta.persistence.*;
import java.time.Instant;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.hibernate.annotations.CreationTimestamp;
import org.hibernate.annotations.SourceType;
import org.hibernate.annotations.UpdateTimestamp;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "user_quota_allocations")
public class UserQuotaAllocation {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_name", nullable = false)
  private PipelinesEnum pipelineName;

  @Column(name = "user_id", nullable = false)
  private String userId;

  @Column(name = "quota_source", nullable = false)
  private QuotaAllocationSourceEnum quotaSource;

  @Column(name = "quota_allocated", nullable = false)
  private int quotaAllocated;

  @Column(name = "quota_consumed", nullable = false)
  private int quotaConsumed;

  @Column(name = "quota_status", nullable = false)
  private QuotaAllocationStatusEnum quotaStatus;

  @Column(name = "created", insertable = false)
  @CreationTimestamp(source = SourceType.DB)
  private Instant created;

  @Column(name = "updated", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private Instant updated;

  @Column(name = "comments")
  private String comments;

  public UserQuotaAllocation(
      PipelinesEnum pipelineName,
      String userId,
      QuotaAllocationSourceEnum quotaSource,
      int quotaAllocated,
      int quotaConsumed,
      QuotaAllocationStatusEnum quotaStatus,
      String comments) {
    this.pipelineName = pipelineName;
    this.userId = userId;
    this.quotaSource = quotaSource;
    this.quotaAllocated = quotaAllocated;
    this.quotaConsumed = quotaConsumed;
    this.quotaStatus = quotaStatus;
    this.comments = comments;
  }
}
