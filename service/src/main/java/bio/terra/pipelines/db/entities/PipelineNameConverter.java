package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class PipelineNameConverter implements AttributeConverter<PipelinesEnum, String> {
  @Override
  public String convertToDatabaseColumn(PipelinesEnum pipelineNameEnum) {
    return pipelineNameEnum.getValue();
  }

  @Override
  public PipelinesEnum convertToEntityAttribute(String pipelineNameString) {
    return PipelinesEnum.valueOf(pipelineNameString.toUpperCase());
  }
}
