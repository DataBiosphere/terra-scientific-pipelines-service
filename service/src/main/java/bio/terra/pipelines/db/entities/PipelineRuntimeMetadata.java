package bio.terra.pipelines.db.entities;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import java.time.Instant;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_runtime_metadata")
public class PipelineRuntimeMetadata {
  @Id
  @Column(name = "pipeline_key", nullable = false)
  private String pipelineKey;

  @Column(name = "pipeline_name")
  private String pipelineName;

  @Column(name = "tool_version")
  private String toolVersion;

  @Column(name = "workspace_billing_project")
  private String workspaceBillingProject;

  @Column(name = "workspace_name")
  private String workspaceName;

  @Column(name = "workspace_google_project")
  private String workspaceGoogleProject;

  @Column(name = "workspace_storage_container_name")
  private String workspaceStorageContainerName;

  @Column(name = "hidden", nullable = false)
  private boolean hidden;

  @Column(name = "updated", insertable = false, updatable = false)
  private Instant updated;
}
