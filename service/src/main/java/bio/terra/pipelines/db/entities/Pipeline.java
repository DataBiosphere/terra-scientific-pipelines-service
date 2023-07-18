package bio.terra.pipelines.db.entities;

import java.util.StringJoiner;
import javax.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

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
}
