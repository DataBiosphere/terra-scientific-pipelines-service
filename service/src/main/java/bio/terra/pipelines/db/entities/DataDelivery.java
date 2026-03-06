package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import java.time.Instant;
import java.util.UUID;
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
@Table(name = "data_delivery")
public class DataDelivery {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_run_id", nullable = false)
  private Long pipelineRunId;

  @ManyToOne(fetch = FetchType.LAZY)
  @JoinColumn(
      name = "pipeline_run_id",
      referencedColumnName = "id",
      insertable = false,
      updatable = false)
  private PipelineRun pipelineRun;

  @Column(name = "job_id", nullable = false)
  private UUID jobId;

  @Column(name = "status", nullable = false)
  private String status;

  @Column(name = "gcs_destination_path", nullable = false)
  private String gcsDestinationPath;

  @Column(name = "created", insertable = false)
  @CreationTimestamp(source = SourceType.DB)
  private Instant created;

  @Column(name = "updated", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private Instant updated;

  public DataDelivery(Long pipelineRunId, UUID jobId, String status, String gcsDestinationPath) {
    this.pipelineRunId = pipelineRunId;
    this.jobId = jobId;
    this.status = status;
    this.gcsDestinationPath = gcsDestinationPath;
  }
}
