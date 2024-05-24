package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import java.util.List;
import java.util.StringJoiner;
import java.util.UUID;
import java.util.stream.Collectors;
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
  private String name;

  @Column(name = "version", nullable = false)
  private String version;

  @Column(name = "display_name", nullable = false)
  private String displayName;

  @Column(name = "description")
  private String description;

  @Column(name = "pipeline_type")
  private String pipelineType;

  @Column(name = "wdl_url")
  private String wdlUrl;

  @Column(name = "wdl_method_name")
  private String wdlMethodName;

  @Column(name = "workspace_id")
  private UUID workspaceId;

  // Note: we fetch eagerly despite not always needing inputs definitions because
  // the number of inputs definitions is expected to be small. Beware using OneToMany with
  // eager fetch on large collections.
  @OneToMany(mappedBy = "pipelineId", fetch = FetchType.EAGER)
  private List<PipelineInputDefinition> pipelineInputDefinitions;

  public Pipeline(
      String name,
      String version,
      String displayName,
      String description,
      String pipelineType,
      String wdlUrl,
      String wdlMethodName,
      UUID workspaceId,
      List<PipelineInputDefinition> pipelineInputDefinitions) {
    this.name = name;
    this.version = version;
    this.displayName = displayName;
    this.description = description;
    this.pipelineType = pipelineType;
    this.wdlUrl = wdlUrl;
    this.wdlMethodName = wdlMethodName;
    this.workspaceId = workspaceId;
    this.pipelineInputDefinitions = pipelineInputDefinitions;
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", Pipeline.class.getSimpleName() + "[", "]")
        .add("pipelineName=" + name)
        .add("version=" + version)
        .add("displayName=" + displayName)
        .add("description=" + description)
        .add("pipelineType=" + pipelineType)
        .add("wdlUrl=" + wdlUrl)
        .add("wdlMethodName=" + wdlMethodName)
        .add("workspaceId=" + workspaceId)
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
        .append(displayName)
        .append(description)
        .append(pipelineType)
        .append(wdlUrl)
        .append(wdlMethodName)
        .append(workspaceId)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof Pipeline)) return false;
    if (obj == this) return true;

    Pipeline otherObject = (Pipeline) obj;
    return new EqualsBuilder()
        .append(id, otherObject.id)
        .append(name, otherObject.name)
        .append(version, otherObject.version)
        .append(displayName, otherObject.displayName)
        .append(description, otherObject.description)
        .append(pipelineType, otherObject.pipelineType)
        .append(wdlUrl, otherObject.wdlUrl)
        .append(wdlMethodName, otherObject.wdlMethodName)
        .append(workspaceId, otherObject.workspaceId)
        .isEquals();
  }

  /**
   * Get a copy of the pipeline input definitions This is required for this object to be stored
   * properly in the working map. See second answer here:
   * https://stackoverflow.com/questions/15833979/java-jackson-deserialize-complex-polymorphic-object-model-jsonmappingexception
   */
  @SuppressWarnings("java:S6204") // we do need to stream then collect
  public List<PipelineInputDefinition> getPipelineInputDefinitions() {
    return pipelineInputDefinitions.stream().collect(Collectors.toList());
  }
}
