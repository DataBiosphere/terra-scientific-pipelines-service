package bio.terra.pipelines.db.entities.converters;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class PipelineVariableTypeConverter
    implements AttributeConverter<PipelineVariableTypesEnum, String> {
  @Override
  public String convertToDatabaseColumn(PipelineVariableTypesEnum pipelineVariableTypeEnum) {
    return pipelineVariableTypeEnum.toString();
  }

  @Override
  public PipelineVariableTypesEnum convertToEntityAttribute(String pipelineVariableTypeString) {
    return PipelineVariableTypesEnum.valueOf(pipelineVariableTypeString);
  }
}
