package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import java.time.LocalDateTime;
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
@Table(name = "pipeline_runs")
@SuppressWarnings({"java:S107"}) // Disable "Methods should not have too many parameters"
public class PipelineRun {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "job_id", nullable = false, unique = true)
  private UUID jobId;

  @Column(name = "user_id", nullable = false)
  private String userId;

  @Column(name = "pipeline_id", nullable = false)
  private Long pipelineId;

  @Column(name = "workspace_id", nullable = false)
  private UUID workspaceId;

  @Column(name = "workspace_storage_container_url", nullable = false)
  private String workspaceStorageContainerUrl;

  @Column(name = "created", insertable = false)
  @CreationTimestamp(source = SourceType.DB)
  private LocalDateTime created;

  @Column(name = "updated", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private LocalDateTime updated;

  @Column(name = "status", nullable = false)
  private CommonPipelineRunStatusEnum status;

  @Column(name = "description")
  private String description;

  @Column(name = "result_url")
  private String resultUrl;

  @Column(name = "is_success")
  private Boolean isSuccess;

  /** Constructor for in progress or complete PipelineRun. */
  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      UUID workspaceId,
      String workspaceStorageContainerUrl,
      LocalDateTime created,
      LocalDateTime updated,
      CommonPipelineRunStatusEnum status,
      String description,
      String resultUrl,
      Boolean isSuccess) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.workspaceId = workspaceId;
    this.workspaceStorageContainerUrl = workspaceStorageContainerUrl;
    this.created = created;
    this.updated = updated;
    this.status = status;
    this.description = description;
    this.resultUrl = resultUrl;
    this.isSuccess = isSuccess;
  }

  /** Constructor for creating a new pipeline run. Timestamps are auto-generated. */
  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      UUID workspaceId,
      String workspaceStorageContainerUrl,
      CommonPipelineRunStatusEnum status) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.workspaceId = workspaceId;
    this.workspaceStorageContainerUrl = workspaceStorageContainerUrl;
    this.status = status;
  }
}
