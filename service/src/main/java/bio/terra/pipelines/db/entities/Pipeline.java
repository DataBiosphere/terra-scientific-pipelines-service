package bio.terra.pipelines.db.entities;

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
    name = "pipelines",
    uniqueConstraints = {@UniqueConstraint(columnNames = {"name", "version"})})
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

  public Pipeline(
      String name,
      String version,
      String displayName,
      String description,
      String pipelineType,
      String wdlUrl,
      String wdlMethodName) {
    this.name = name;
    this.version = version;
    this.displayName = displayName;
    this.description = description;
    this.pipelineType = pipelineType;
    this.wdlUrl = wdlUrl;
    this.wdlMethodName = wdlMethodName;
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
        .isEquals();
  }
}
