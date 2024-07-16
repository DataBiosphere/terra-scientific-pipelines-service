package bio.terra.pipelines.db.entities.converters;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class PipelineInputTypeConverter
    implements AttributeConverter<PipelineVariableTypesEnum, String> {
  @Override
  public String convertToDatabaseColumn(PipelineVariableTypesEnum pipelineInputTypeEnum) {
    return pipelineInputTypeEnum.toString();
  }

  @Override
  public PipelineVariableTypesEnum convertToEntityAttribute(String pipelineInputTypeString) {
    return PipelineVariableTypesEnum.valueOf(pipelineInputTypeString);
  }
}
