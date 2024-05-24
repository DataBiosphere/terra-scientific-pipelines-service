package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_input_definitions")
public class PipelineInputDefinition {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_id")
  private Long pipelineId;

  @Column(name = "name", nullable = false)
  private String name;

  @Column(name = "type", nullable = false)
  private PipelineInputTypesEnum type;

  @Column(name = "is_required", nullable = false)
  private Boolean isRequired;

  @Column(name = "user_provided", nullable = false)
  private Boolean userProvided;

  @Column(name = "default_value")
  private String defaultValue; // must be a String representation of the value

  @JsonCreator
  public PipelineInputDefinition(
      @JsonProperty("pipelineId") Long pipelineId,
      @JsonProperty("name") String name,
      @JsonProperty("type") PipelineInputTypesEnum type,
      @JsonProperty("isRequired") Boolean isRequired,
      @JsonProperty("userProvided") Boolean userProvided,
      @JsonProperty("defaultValue") String defaultValue) {
    this.pipelineId = pipelineId;
    this.name = name;
    this.type = type;
    this.isRequired = isRequired;
    this.userProvided = userProvided;
    this.defaultValue = defaultValue;
  }
}
