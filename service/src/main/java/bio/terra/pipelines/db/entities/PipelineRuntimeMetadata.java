package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import jakarta.persistence.*;
import java.util.StringJoiner;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(
    name = "pipeline_runtime_metadata",
    uniqueConstraints = {@UniqueConstraint(columnNames = {"name", "version"})})
@SuppressWarnings("java:S107") // Disable "Methods should not have too many parameters"
public class PipelineRuntimeMetadata {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "name", nullable = false)
  private PipelinesEnum name;

  @Column(name = "version", nullable = false)
  private Integer version;

  @Column(name = "hidden", nullable = false)
  private boolean hidden;

  @Column(name = "tool_version")
  private String toolVersion;

  @Column(name = "workspace_billing_project")
  private String workspaceBillingProject;

  @Column(name = "workspace_name")
  private String workspaceName;

  @Column(name = "workspace_storage_container_name")
  private String workspaceStorageContainerName;

  @Column(name = "workspace_google_project")
  private String workspaceGoogleProject;

  public PipelineRuntimeMetadata(
      PipelinesEnum name,
      Integer version,
      boolean hidden,
      String toolVersion,
      String workspaceBillingProject,
      String workspaceName,
      String workspaceStorageContainerName,
      String workspaceGoogleProject) {
    this.name = name;
    this.version = version;
    this.hidden = hidden;
    this.toolVersion = toolVersion;
    this.workspaceBillingProject = workspaceBillingProject;
    this.workspaceName = workspaceName;
    this.workspaceStorageContainerName = workspaceStorageContainerName;
    this.workspaceGoogleProject = workspaceGoogleProject;
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", PipelineRuntimeMetadata.class.getSimpleName() + "[", "]")
        .add("pipelineName=" + name)
        .add("version=" + version)
        .add("hidden=" + hidden)
        .add("toolVersion=" + toolVersion)
        .add("workspaceBillingProject=" + workspaceBillingProject)
        .add("workspaceName=" + workspaceName)
        .add("workspaceStorageContainerName=" + workspaceStorageContainerName)
        .add("workspaceGoogleProject=" + workspaceGoogleProject)
        .toString();
  }

  @SuppressWarnings("java:S125") // The comment here isn't "commented code"
  // we override equals() below so that we can compare Pipeline objects in tests;
  // according to
  // https://stackoverflow.com/questions/27581/what-issues-should-be-considered-when-overriding-equals-and-hashcode-in-java/27609#27609
  // we should override hashCode() if we override equals()
  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        // two randomly chosen prime numbers
        // if deriving: appendSuper(super.hashCode()).
        .append(id)
        .append(name)
        .append(version)
        .append(hidden)
        .append(toolVersion)
        .append(workspaceBillingProject)
        .append(workspaceName)
        .append(workspaceStorageContainerName)
        .append(workspaceGoogleProject)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PipelineRuntimeMetadata otherObject)) return false;
    if (obj == this) return true;

    return new EqualsBuilder()
        .append(id, otherObject.id)
        .append(name, otherObject.name)
        .append(version, otherObject.version)
        .append(hidden, otherObject.hidden)
        .append(toolVersion, otherObject.toolVersion)
        .append(workspaceBillingProject, otherObject.workspaceBillingProject)
        .append(workspaceName, otherObject.workspaceName)
        .append(workspaceStorageContainerName, otherObject.workspaceStorageContainerName)
        .append(workspaceGoogleProject, otherObject.workspaceGoogleProject)
        .isEquals();
  }
}
