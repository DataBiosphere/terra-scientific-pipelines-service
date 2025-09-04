package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
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
@Table(name = "pipeline_output_definitions")
public class PipelineOutputDefinition extends BasePipelineVariableDefinition {

  public PipelineOutputDefinition(
      Long pipelineId,
      String name,
      String wdlVariableName,
      PipelineVariableTypesEnum type,
      boolean isRequired) {
    super(pipelineId, name, wdlVariableName, type, isRequired);
  }

  @SuppressWarnings("java:S125") // The comment here isn't "commented code"
  // we override equals() below so that we can compare PipelineOutputDefinition objects in tests;
  // according to
  // https://stackoverflow.com/questions/27581/what-issues-should-be-considered-when-overriding-equals-and-hashcode-in-java/27609#27609
  // we should override hashCode() if we override equals()
  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        // two randomly chosen prime numbers
        .append(getId())
        .append(getPipelineId())
        .append(getName())
        .append(getWdlVariableName())
        .append(getType())
        .append(isRequired())
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PipelineOutputDefinition otherObject)) return false;
    if (obj == this) return true;
    return new EqualsBuilder()
        .append(getId(), otherObject.getId())
        .append(getPipelineId(), otherObject.getPipelineId())
        .append(getName(), otherObject.getName())
        .append(getWdlVariableName(), otherObject.getWdlVariableName())
        .append(getType(), otherObject.getType())
        .append(isRequired(), otherObject.isRequired())
        .isEquals();
  }
}
