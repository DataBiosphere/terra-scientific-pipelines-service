package bio.terra.pipelines.db.entities.converters;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class PipelineInputTypeConverter
    implements AttributeConverter<PipelineInputTypesEnum, String> {
  @Override
  public String convertToDatabaseColumn(PipelineInputTypesEnum pipelineInputTypeEnum) {
    return pipelineInputTypeEnum.toString();
  }

  @Override
  public PipelineInputTypesEnum convertToEntityAttribute(String pipelineInputTypeString) {
    return PipelineInputTypesEnum.valueOf(pipelineInputTypeString);
  }
}
