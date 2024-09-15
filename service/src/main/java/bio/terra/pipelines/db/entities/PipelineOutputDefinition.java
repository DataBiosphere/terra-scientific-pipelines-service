package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import jakarta.persistence.Entity;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.apache.commons.lang3.builder.EqualsBuilder;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "pipeline_output_definitions")
public class PipelineOutputDefinition extends BasePipelineVariableDefinition {

  public PipelineOutputDefinition(
      Long pipelineId, String name, String wdlVariableName, PipelineVariableTypesEnum type) {
    super(pipelineId, name, wdlVariableName, type);
  }

  // we override equals() below so that we can compare PipelineInputDefinition objects;
  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PipelineOutputDefinition otherObject)) return false;
    if (obj == this) return true;
    return new EqualsBuilder()
        .append(getPipelineId(), otherObject.getPipelineId())
        .append(getName(), otherObject.getName())
        .append(getWdlVariableName(), otherObject.getWdlVariableName())
        .append(getType(), otherObject.getType())
        .isEquals();
  }
}
