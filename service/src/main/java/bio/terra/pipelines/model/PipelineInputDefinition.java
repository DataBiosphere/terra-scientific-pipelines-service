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
@SuppressWarnings("java:S107")
public class PipelineInputDefinition extends PipelineVariableDefinition {
  private boolean userProvided;
  private String defaultValue;
  private Double minValue;
  private Double maxValue;
  private String fileSuffix;

  public PipelineInputDefinition(
      String pipelineKey,
      String name,
      String wdlVariableName,
      String displayName,
      String description,
      PipelineVariableTypesEnum type,
      String fileSuffix,
      boolean isRequired,
      boolean userProvided,
      String defaultValue,
      Double minValue,
      Double maxValue) {
    super(pipelineKey, name, wdlVariableName, displayName, description, type, isRequired);
    this.fileSuffix = fileSuffix;
    this.userProvided = userProvided;
    this.defaultValue = defaultValue;
    this.minValue = minValue;
    this.maxValue = maxValue;
  }

  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        .appendSuper(super.hashCode())
        .append(fileSuffix)
        .append(userProvided)
        .append(defaultValue)
        .append(minValue)
        .append(maxValue)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PipelineInputDefinition otherObject)) return false;
    if (obj == this) return true;
    return new EqualsBuilder()
        .appendSuper(super.equals(otherObject))
        .append(fileSuffix, otherObject.fileSuffix)
        .append(userProvided, otherObject.userProvided)
        .append(defaultValue, otherObject.defaultValue)
        .append(minValue, otherObject.minValue)
        .append(maxValue, otherObject.maxValue)
        .isEquals();
  }
}
