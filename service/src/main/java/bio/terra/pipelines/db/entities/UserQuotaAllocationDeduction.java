package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import java.time.Instant;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.hibernate.annotations.CreationTimestamp;
import org.hibernate.annotations.SourceType;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "user_quota_allocation_deductions")
public class UserQuotaAllocationDeduction {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "user_quota_allocation_id", nullable = false)
  private Long userQuotaAllocationId;

  @Column(name = "pipeline_run_id")
  private Long pipelineRunId;

  @Column(name = "amount", nullable = false)
  private int amount;

  @Column(name = "comments")
  private String comments;

  @Column(name = "created", insertable = false)
  @CreationTimestamp(source = SourceType.DB)
  private Instant created;

  public UserQuotaAllocationDeduction(
      Long userQuotaAllocationId, Long pipelineRunId, int amount, String comments) {
    this.userQuotaAllocationId = userQuotaAllocationId;
    this.pipelineRunId = pipelineRunId;
    this.amount = amount;
    this.comments = comments;
  }
}
