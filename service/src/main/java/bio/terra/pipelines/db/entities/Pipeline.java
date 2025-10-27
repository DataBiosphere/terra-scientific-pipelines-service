package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import jakarta.persistence.*;
import java.util.ArrayList;
import java.util.List;
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
    name = "pipelines",
    uniqueConstraints = {@UniqueConstraint(columnNames = {"name", "version"})})
@SuppressWarnings("java:S107") // Disable "Methods should not have too many parameters"
public class Pipeline {
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

  @Column(name = "display_name", nullable = false)
  private String displayName;

  @Column(name = "description")
  private String description;

  @Column(name = "pipeline_type")
  private String pipelineType;

  @Column(name = "wdl_url")
  private String wdlUrl;

  @Column(name = "tool_name")
  private String toolName;

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

  // Note: we fetch eagerly despite not always needing inputs definitions because
  // the number of inputs definitions is expected to be small. Beware using OneToMany with
  // eager fetch on large collections.
  @OneToMany(mappedBy = "pipelineId", fetch = FetchType.EAGER)
  private List<PipelineInputDefinition> pipelineInputDefinitions;

  @OneToMany(mappedBy = "pipelineId", fetch = FetchType.EAGER)
  private List<PipelineOutputDefinition> pipelineOutputDefinitions;

  public Pipeline(
      PipelinesEnum name,
      Integer version,
      boolean hidden,
      String displayName,
      String description,
      String pipelineType,
      String wdlUrl,
      String toolName,
      String toolVersion,
      String workspaceBillingProject,
      String workspaceName,
      String workspaceStorageContainerName,
      String workspaceGoogleProject,
      List<PipelineInputDefinition> pipelineInputDefinitions,
      List<PipelineOutputDefinition> pipelineOutputDefinitions) {
    this.name = name;
    this.version = version;
    this.hidden = hidden;
    this.displayName = displayName;
    this.description = description;
    this.pipelineType = pipelineType;
    this.wdlUrl = wdlUrl;
    this.toolName = toolName;
    this.toolVersion = toolVersion;
    this.workspaceBillingProject = workspaceBillingProject;
    this.workspaceName = workspaceName;
    this.workspaceStorageContainerName = workspaceStorageContainerName;
    this.workspaceGoogleProject = workspaceGoogleProject;
    this.pipelineInputDefinitions = pipelineInputDefinitions;
    this.pipelineOutputDefinitions = pipelineOutputDefinitions;
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", Pipeline.class.getSimpleName() + "[", "]")
        .add("pipelineName=" + name)
        .add("version=" + version)
        .add("hidden=" + hidden)
        .add("displayName=" + displayName)
        .add("description=" + description)
        .add("pipelineType=" + pipelineType)
        .add("wdlUrl=" + wdlUrl)
        .add("toolName=" + toolName)
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
        .append(displayName)
        .append(description)
        .append(pipelineType)
        .append(wdlUrl)
        .append(toolName)
        .append(toolVersion)
        .append(workspaceBillingProject)
        .append(workspaceName)
        .append(workspaceStorageContainerName)
        .append(workspaceGoogleProject)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof Pipeline otherObject)) return false;
    if (obj == this) return true;

    return new EqualsBuilder()
        .append(id, otherObject.id)
        .append(name, otherObject.name)
        .append(version, otherObject.version)
        .append(hidden, otherObject.hidden)
        .append(displayName, otherObject.displayName)
        .append(description, otherObject.description)
        .append(pipelineType, otherObject.pipelineType)
        .append(wdlUrl, otherObject.wdlUrl)
        .append(toolName, otherObject.toolName)
        .append(toolVersion, otherObject.toolVersion)
        .append(workspaceBillingProject, otherObject.workspaceBillingProject)
        .append(workspaceName, otherObject.workspaceName)
        .append(workspaceStorageContainerName, otherObject.workspaceStorageContainerName)
        .append(workspaceGoogleProject, otherObject.workspaceGoogleProject)
        .isEquals();
  }

  /**
   * Get a copy of the pipeline input definitions This is required for this object to be properly
   * deserialized from the Stairway working map. See second answer here:
   * https://stackoverflow.com/questions/15833979/java-jackson-deserialize-complex-polymorphic-object-model-jsonmappingexception
   */
  public List<PipelineInputDefinition> getPipelineInputDefinitions() {
    return new ArrayList<>(pipelineInputDefinitions);
  }

  /**
   * Get a copy of the pipeline output definitions, in a similar fashion to
   * getPipelineInputDefinitions().
   */
  public List<PipelineOutputDefinition> getPipelineOutputDefinitions() {
    return new ArrayList<>(pipelineOutputDefinitions);
  }
}
