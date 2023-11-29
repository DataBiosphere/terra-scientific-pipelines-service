package bio.terra.pipelines.db.entities;

import java.util.StringJoiner;
import javax.persistence.*;
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
    uniqueConstraints = {@UniqueConstraint(columnNames = {"pipeline_id", "version"})})
public class Pipeline {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_id", nullable = false)
  private String pipelineId;

  @Column(name = "version", nullable = false)
  private String version;

  @Column(name = "display_name", nullable = false)
  private String displayName;

  @Column(name = "description")
  private String description;

  public Pipeline(String pipelineId, String version, String displayName, String description) {
    this.pipelineId = pipelineId;
    this.version = version;
    this.displayName = displayName;
    this.description = description;
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", Pipeline.class.getSimpleName() + "[", "]")
        .add("pipelineId=" + pipelineId)
        .add("version=" + version)
        .add("displayName=" + displayName)
        .add("description=" + description)
        .toString();
  }

  // we override equals() below so that we can compare Pipeline objects in tests;
  // according to https://stackoverflow.com/questions/27581/what-issues-should-be-considered-when-overriding-equals-and-hashcode-in-java/27609#27609
  // we should override hashCode() if we override equals()
  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        . // two randomly chosen prime numbers
        // if deriving: appendSuper(super.hashCode()).
        append(pipelineId)
        .append(version)
        .append(displayName)
        .append(description)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof Pipeline)) return false;
    if (obj == this) return true;

    Pipeline otherObject = (Pipeline) obj;
    return new EqualsBuilder()
        .append(pipelineId, otherObject.pipelineId)
        .append(version, otherObject.version)
        .append(displayName, otherObject.displayName)
        .append(description, otherObject.description)
        .isEquals();
  }
}
