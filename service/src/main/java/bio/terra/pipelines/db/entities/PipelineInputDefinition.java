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

  @Column(name = "user_provided", nullable = false)
  private boolean userProvided;

  @Column(name = "expects_custom_value", nullable = false)
  private boolean expectsCustomValue; // only relevant for service-provided inputs

  @Column(name = "default_value")
  private String defaultValue; // must be a String representation of the value

  @Column(name = "min_value")
  private Integer minValue; // inclusive minimum

  @Column(name = "max_value")
  private Integer maxValue; // inclusive maximum

  @Column(name = "file_suffix")
  private String fileSuffix;

  public PipelineInputDefinition(
      Long pipelineId,
      String name,
      String wdlVariableName,
      PipelineVariableTypesEnum type,
      String fileSuffix,
      boolean isRequired,
      boolean userProvided,
      boolean expectsCustomValue,
      String defaultValue,
      Integer minValue,
      Integer maxValue) {
    super(pipelineId, name, wdlVariableName, type, isRequired);
    this.fileSuffix = fileSuffix;
    this.userProvided = userProvided;
    this.expectsCustomValue = expectsCustomValue;
    this.defaultValue = defaultValue;
    this.minValue = minValue;
    this.maxValue = maxValue;
  }

  @SuppressWarnings("java:S125") // The comment here isn't "commented code"
  // we override equals() below so that we can compare PipelineInputDefinition objects in tests;
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
        .append(fileSuffix)
        .append(isRequired())
        .append(userProvided)
        .append(expectsCustomValue)
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
        .append(getId(), otherObject.getId())
        .append(getPipelineId(), otherObject.getPipelineId())
        .append(getName(), otherObject.getName())
        .append(getWdlVariableName(), otherObject.getWdlVariableName())
        .append(getType(), otherObject.getType())
        .append(fileSuffix, otherObject.fileSuffix)
        .append(isRequired(), otherObject.isRequired())
        .append(userProvided, otherObject.userProvided)
        .append(expectsCustomValue, otherObject.expectsCustomValue)
        .append(defaultValue, otherObject.defaultValue)
        .append(minValue, otherObject.minValue)
        .append(maxValue, otherObject.maxValue)
        .isEquals();
  }
}
