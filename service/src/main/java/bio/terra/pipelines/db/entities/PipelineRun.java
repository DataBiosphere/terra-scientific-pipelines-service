package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
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

  @ManyToOne(fetch = FetchType.LAZY)
  @JoinColumn(
      name = "pipeline_id",
      referencedColumnName = "id",
      insertable = false,
      updatable = false)
  private Pipeline pipeline;

  @Column(name = "tool_version")
  private String toolVersion;

  @Column(name = "workspace_billing_project")
  private String workspaceBillingProject;

  @Column(name = "workspace_name")
  private String workspaceName;

  @Column(name = "workspace_storage_container_name", nullable = false)
  private String workspaceStorageContainerName;

  @Column(name = "workspace_google_project", nullable = false)
  private String workspaceGoogleProject;

  @Column(name = "created", insertable = false)
  @CreationTimestamp(source = SourceType.DB)
  private Instant created;

  @Column(name = "updated", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private Instant updated;

  @Column(name = "status", nullable = false)
  private CommonPipelineRunStatusEnum status;

  @Column(name = "description")
  private String description;

  @Column(name = "quota_consumed")
  private Integer quotaConsumed;

  @Column(name = "raw_quota_consumed")
  private Integer rawQuotaConsumed;

  /** Constructor for in progress or complete PipelineRun. */
  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      String toolVersion,
      String workspaceBillingProject,
      String workspaceName,
      String workspaceStorageContainerName,
      String workspaceGoogleProject,
      Instant created,
      Instant updated,
      CommonPipelineRunStatusEnum status,
      String description,
      Integer quotaConsumed,
      Integer rawQuotaConsumed) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.toolVersion = toolVersion;
    this.workspaceBillingProject = workspaceBillingProject;
    this.workspaceName = workspaceName;
    this.workspaceStorageContainerName = workspaceStorageContainerName;
    this.workspaceGoogleProject = workspaceGoogleProject;
    this.created = created;
    this.updated = updated;
    this.status = status;
    this.description = description;
    this.quotaConsumed = quotaConsumed;
    this.rawQuotaConsumed = rawQuotaConsumed;
  }

  /** Constructor for creating a new pipeline run. Timestamps are auto-generated. */
  public PipelineRun(
      UUID jobId,
      String userId,
      Long pipelineId,
      String toolVersion,
      String workspaceBillingProject,
      String workspaceName,
      String workspaceStorageContainerName,
      String workspaceGoogleProject,
      CommonPipelineRunStatusEnum status,
      String description) {
    this.jobId = jobId;
    this.userId = userId;
    this.pipelineId = pipelineId;
    this.toolVersion = toolVersion;
    this.workspaceBillingProject = workspaceBillingProject;
    this.workspaceName = workspaceName;
    this.workspaceStorageContainerName = workspaceStorageContainerName;
    this.workspaceGoogleProject = workspaceGoogleProject;
    this.status = status;
    this.description = description;
  }
}
