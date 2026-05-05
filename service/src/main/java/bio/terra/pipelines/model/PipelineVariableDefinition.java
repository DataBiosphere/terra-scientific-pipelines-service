package bio.terra.pipelines.model;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

@Getter
@Setter
@NoArgsConstructor
public abstract class PipelineVariableDefinition {
  private String pipelineKey;
  private String name;
  private String wdlVariableName;
  private String displayName;
  private String description;
  private PipelineVariableTypesEnum type;
  private boolean required;

  protected PipelineVariableDefinition(
      String pipelineKey,
      String name,
      String wdlVariableName,
      String displayName,
      String description,
      PipelineVariableTypesEnum type,
      boolean required) {
    this.pipelineKey = pipelineKey;
    this.name = name;
    this.wdlVariableName = wdlVariableName;
    this.displayName = displayName;
    this.description = description;
    this.type = type;
    this.required = required;
  }

  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        .append(pipelineKey)
        .append(name)
        .append(wdlVariableName)
        .append(displayName)
        .append(description)
        .append(type)
        .append(required)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PipelineVariableDefinition otherObject)) return false;
    if (obj == this) return true;
    return new EqualsBuilder()
        .append(pipelineKey, otherObject.pipelineKey)
        .append(name, otherObject.name)
        .append(wdlVariableName, otherObject.wdlVariableName)
        .append(displayName, otherObject.displayName)
        .append(description, otherObject.description)
        .append(type, otherObject.type)
        .append(required, otherObject.required)
        .isEquals();
  }
}
