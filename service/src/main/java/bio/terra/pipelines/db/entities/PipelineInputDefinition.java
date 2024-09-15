package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_input_definitions")
@SuppressWarnings("java:S107") // Disable "Methods should not have too many parameters"
public class PipelineInputDefinition extends BasePipelineVariableDefinition {
  @Column(name = "is_required", nullable = false)
  private Boolean isRequired;

  @Column(name = "user_provided", nullable = false)
  private Boolean userProvided;

  @Column(name = "default_value")
  private String defaultValue; // must be a String representation of the value

  @Column(name = "file_suffix")
  private String fileSuffix;

  public PipelineInputDefinition(
      Long pipelineId,
      String name,
      String wdlVariableName,
      PipelineVariableTypesEnum type,
      String fileSuffix,
      Boolean isRequired,
      Boolean userProvided,
      String defaultValue) {
    super(pipelineId, name, wdlVariableName, type);
    this.fileSuffix = fileSuffix;
    this.isRequired = isRequired;
    this.userProvided = userProvided;
    this.defaultValue = defaultValue;
  }

  // we override equals() below so that we can compare PipelineInputDefinition objects;
  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        // two randomly chosen prime numbers
        .append(fileSuffix)
        .append(isRequired)
        .append(userProvided)
        .append(defaultValue)
        .appendSuper(super.hashCode())
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PipelineInputDefinition otherObject)) return false;
    if (obj == this) return true;
    return new EqualsBuilder()
        .append(fileSuffix, otherObject.fileSuffix)
        .append(isRequired, otherObject.isRequired)
        .append(userProvided, otherObject.userProvided)
        .append(defaultValue, otherObject.defaultValue)
        .append(getPipelineId(), otherObject.getPipelineId())
        .append(getName(), otherObject.getName())
        .append(getWdlVariableName(), otherObject.getWdlVariableName())
        .append(getType(), otherObject.getType())
        .isEquals();
  }
}
