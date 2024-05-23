package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.PipelineInputTypesEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;
import java.util.stream.Stream;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class PipelineInputTypeConverter
    implements AttributeConverter<PipelineInputTypesEnum, String> {
  @Override
  public String convertToDatabaseColumn(PipelineInputTypesEnum pipelineInputTypeEnum) {
    if (pipelineInputTypeEnum == null) {
      return null;
    }
    return pipelineInputTypeEnum.toString();
  }

  @Override
  public PipelineInputTypesEnum convertToEntityAttribute(String pipelineInputTypeString) {
    if (pipelineInputTypeString == null) {
      return null;
    }

    return Stream.of(PipelineInputTypesEnum.values())
        .filter(c -> c.toString().equals(pipelineInputTypeString))
        .findFirst()
        .orElseThrow(IllegalArgumentException::new);
  }
}
