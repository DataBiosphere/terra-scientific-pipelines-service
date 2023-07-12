package bio.terra.pipelines.db.entities;

import java.util.StringJoiner;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipelines")
public class Pipeline {
  @Id
  @Column(name = "pipeline_id", nullable = false)
  private String pipelineId;

  @Column(name = "display_name")
  private String displayName;

  @Column(name = "description")
  private String description;

  public Pipeline(String pipelineId, String displayName, String description) {
    this.pipelineId = pipelineId;
    this.displayName = displayName;
    this.description = description;
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", Pipeline.class.getSimpleName() + "[", "]")
        .add("pipelineId=" + pipelineId)
        .add("displayName=" + displayName)
        .add("description=" + description)
        .toString();
  }
}
