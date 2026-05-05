package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelineKeyUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

@Getter
@Setter
@NoArgsConstructor
@SuppressWarnings("java:S107")
public class Pipeline {
  private PipelinesEnum name;
  private Integer version;
  private boolean hidden;
  private String displayName;
  private String description;
  private String pipelineType;
  private String toolName;
  private String toolVersion;
  private String workspaceBillingProject;
  private String workspaceName;
  private String workspaceStorageContainerName;
  private String workspaceGoogleProject;
  private List<PipelineInputDefinition> pipelineInputDefinitions;
  private List<PipelineOutputDefinition> pipelineOutputDefinitions;

  public Pipeline(
      PipelinesEnum name,
      Integer version,
      boolean hidden,
      String displayName,
      String description,
      String pipelineType,
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
    this.toolName = toolName;
    this.toolVersion = toolVersion;
    this.workspaceBillingProject = workspaceBillingProject;
    this.workspaceName = workspaceName;
    this.workspaceStorageContainerName = workspaceStorageContainerName;
    this.workspaceGoogleProject = workspaceGoogleProject;
    this.pipelineInputDefinitions = pipelineInputDefinitions;
    this.pipelineOutputDefinitions = pipelineOutputDefinitions;
  }

  public Pipeline(
      PipelinesEnum name,
      Integer version,
      String displayName,
      String description,
      String pipelineType,
      String toolName,
      List<PipelineInputDefinition> pipelineInputDefinitions,
      List<PipelineOutputDefinition> pipelineOutputDefinitions) {
    this.name = name;
    this.version = version;
    this.displayName = displayName;
    this.description = description;
    this.pipelineType = pipelineType;
    this.toolName = toolName;
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
        .add("toolName=" + toolName)
        .add("toolVersion=" + toolVersion)
        .add("workspaceBillingProject=" + workspaceBillingProject)
        .add("workspaceName=" + workspaceName)
        .add("workspaceStorageContainerName=" + workspaceStorageContainerName)
        .add("workspaceGoogleProject=" + workspaceGoogleProject)
        .toString();
  }

  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        .append(name)
        .append(version)
        .append(hidden)
        .append(displayName)
        .append(description)
        .append(pipelineType)
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
        .append(name, otherObject.name)
        .append(version, otherObject.version)
        .append(hidden, otherObject.hidden)
        .append(displayName, otherObject.displayName)
        .append(description, otherObject.description)
        .append(pipelineType, otherObject.pipelineType)
        .append(toolName, otherObject.toolName)
        .append(toolVersion, otherObject.toolVersion)
        .append(workspaceBillingProject, otherObject.workspaceBillingProject)
        .append(workspaceName, otherObject.workspaceName)
        .append(workspaceStorageContainerName, otherObject.workspaceStorageContainerName)
        .append(workspaceGoogleProject, otherObject.workspaceGoogleProject)
        .isEquals();
  }

  /**
   * Get a copy of the pipeline input definitions. This is required for this object to be properly
   * deserialized from the Stairway working map.
   */
  public List<PipelineInputDefinition> getPipelineInputDefinitions() {
    return pipelineInputDefinitions == null
        ? new ArrayList<>()
        : new ArrayList<>(pipelineInputDefinitions);
  }

  /** Get a copy of the pipeline output definitions. */
  public List<PipelineOutputDefinition> getPipelineOutputDefinitions() {
    return pipelineOutputDefinitions == null
        ? new ArrayList<>()
        : new ArrayList<>(pipelineOutputDefinitions);
  }

  /** Canonical pipeline key in the form {pipelineName}_v{pipelineVersion}. */
  public String getKey() {
    return PipelineKeyUtils.buildPipelineKey(name, version);
  }
}
