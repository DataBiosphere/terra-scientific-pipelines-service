package bio.terra.pipelines.db.entities;

import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class PipelineRunStatusConverter
    implements AttributeConverter<CommonPipelineRunStatusEnum, String> {
  @Override
  public String convertToDatabaseColumn(CommonPipelineRunStatusEnum pipelineNameEnum) {
    return pipelineNameEnum.toString();
  }

  @Override
  public CommonPipelineRunStatusEnum convertToEntityAttribute(String pipelineRunStatusString) {
    return CommonPipelineRunStatusEnum.valueOf(pipelineRunStatusString);
  }
}
